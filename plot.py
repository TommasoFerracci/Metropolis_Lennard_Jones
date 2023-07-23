import matplotlib.pyplot as plt 
import pandas as pd

plt.style.use("default")
plt.rc("text", usetex=True)
plt.rc("font", family="cm")
plt.rcParams["grid.color"] = (0.5, 0.5, 0.5, 0.2)

correlation = pd.read_csv('correlation.txt', sep='\t', header=None)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6), dpi=200)
ax.set_xlim([0, 6.5])
ax.set_ylim([0, 3.5])
ax.set_xlabel("r/$\sigma$")
ax.set_ylabel("g(r)")
ax.grid()
ax.plot(correlation[0], correlation[1], linestyle="-", color="purple", linewidth=1)
fig.savefig("images/pair_correlation.png", dpi=200, bbox_inches="tight", pad_inches=0.3)