import os
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("plots", exist_ok=True)

# 1) Insert size distribution
sizes = np.random.lognormal(mean=4.5, sigma=0.35, size=3000).astype(int)
sizes = sizes[(sizes > 50) & (sizes < 800)]
plt.figure(figsize=(8,5))
plt.hist(sizes, bins=50, color="#5DADE2")
plt.title("Insert size distribution")
plt.xlabel("Fragment length (bp)")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("plots/insert_size_distribution.png")
plt.close()

# 2) Peak length distribution
lengths = np.random.lognormal(mean=5.4, sigma=0.55, size=2000).astype(int)
plt.figure(figsize=(8,5))
plt.hist(lengths, bins=50, color="#7DCEA0")
plt.title("Peak length distribution")
plt.xlabel("Length (bp)")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("plots/peak_length_distribution.png")
plt.close()

print("Saved plots: insert_size_distribution.png, peak_length_distribution.png")
