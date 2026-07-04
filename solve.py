import numpy as np
import pandas as pd
from scipy.optimize import least_squares

DATA = "xy_data.csv"

# bounds from the problem statement: theta in (0, 50deg), M in (-0.05, 0.05), X in (0, 100)
LO = np.array([0.0, -0.05, 0.0])
HI = np.array([np.deg2rad(50.0), 0.05, 100.0])

df = pd.read_csv(DATA)
x = df["x"].to_numpy()
y = df["y"].to_numpy()


def residuals(p, x, y):
    # undo the rotation/shift; if p is right, v matches e^(Mt) sin(0.3t)
    theta, M, X = p
    c, s = np.cos(theta), np.sin(theta)
    dx, dy = x - X, y - 42.0
    t = dx * c + dy * s
    v = -dx * s + dy * c
    return v - np.exp(M * np.abs(t)) * np.sin(0.3 * t)


# sin(0.3t) has period ~21, so multi-start or the fit can settle a period off
starts = []
for th_deg in np.arange(2.0, 50.0, 3.0):
    th = np.deg2rad(th_deg)
    # x cos + (y-42) sin = t + X cos(th), and t should span roughly (6, 60)
    u0 = x * np.cos(th) + (y - 42.0) * np.sin(th)
    X0 = ((u0.min() - 6.0) + (u0.max() - 60.0)) / 2.0 / np.cos(th)
    if 0.0 < X0 < 100.0:
        starts.append((th, 0.0, X0))
    for X0 in np.arange(5.0, 100.0, 10.0):
        starts.append((th, 0.0, X0))

best = None
for p0 in starts:
    sol = least_squares(residuals, p0, args=(x, y), bounds=(LO, HI))
    if best is None or sol.cost < best.cost:
        best = sol

theta, M, X = best.x
r = residuals(best.x, x, y)

c, s = np.cos(theta), np.sin(theta)
t = (x - X) * c + (y - 42.0) * s

print(f"theta = {theta:.10f} rad = {np.degrees(theta):.6f} deg")
print(f"M     = {M:.10f}")
print(f"X     = {X:.10f}")
print(f"residual rms {np.sqrt(np.mean(r ** 2)):.3e}, max {np.abs(r).max():.3e}")
print(f"recovered t range [{t.min():.4f}, {t.max():.4f}], expected inside (6, 60)")

# rebuild the curve at the recovered t values and compare with the data
w = np.exp(M * np.abs(t)) * np.sin(0.3 * t)
xf = t * c - w * s + X
yf = 42.0 + t * s + w * c
print(f"L1 over all 1500 points: {np.abs(xf - x).sum() + np.abs(yf - y).sum():.3e}")

print("\nDesmos / latex:")
print(f"\\left(t\\cos({theta:.6f})-e^{{{M:.6f}\\left|t\\right|}}"
      f"\\sin(0.3t)\\sin({theta:.6f})+{X:.6f},\\ "
      f"42+t\\sin({theta:.6f})+e^{{{M:.6f}\\left|t\\right|}}"
      f"\\sin(0.3t)\\cos({theta:.6f})\\right)")
