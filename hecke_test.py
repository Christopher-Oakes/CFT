import cmath
import math
import matplotlib.pyplot as plt

# ==============================================================
# 1. Data: paste from your output (p, q, u, mod)
# ==============================================================

raw_data = [
    (73, 103, 5483, 7519),
    (73, 107, 4476, 7811),
    (73, 199, 7626, 14527),
    (73, 317, 12182, 23141),
    (73, 421, 291, 30733),
    (73, 509, 17635, 37157),
    (73, 521, 25984, 38033),
    (103, 107, 5585, 11021),
    (103, 199, 7651, 20497),
    (103, 317, 19428, 32651),
    (103, 421, 4795, 43363),
    (103, 509, 43920, 52427),
    (103, 521, 48864, 53663),
    (107, 199, 7288, 21293),
    (107, 317, 6972, 33919),
    (107, 421, 1517, 45047),
    (107, 509, 11794, 54463),
    (107, 521, 49607, 55747),
    (199, 317, 29161, 63083),
    (199, 421, 35806, 83779),
    (199, 509, 14412, 101291),
    (199, 521, 89223, 103679),
    (317, 421, 78727, 133457),
    (317, 509, 43916, 161353),
    (317, 521, 11744, 165157),
    (421, 509, 74082, 214289),
    (421, 521, 188324, 219341),
    (509, 521, 134751, 265189)
]

# ==============================================================
# 2. Compute phase embeddings
# ==============================================================

def phase_embedding(u, mod):
    """Map residue u mod m to a complex phase e^{2πi*u/m}."""
    return cmath.exp(2j * math.pi * (u % mod) / mod)

def hecke_phase_deviation(p, q, u, mod):
    """Compute phase deviation between composite and product phases."""
    # Base phases
    phi_p = phase_embedding(p, mod)
    phi_q = phase_embedding(q, mod)
    phi_pq = phase_embedding(u, mod)
    # Phase deviation (radians, principal value)
    delta = abs(cmath.phase(phi_pq / (phi_p * phi_q)))
    return delta

# ==============================================================
# 3. Analyze all pairs
# ==============================================================

deviations = []
for (p, q, u, mod) in raw_data:
    delta = hecke_phase_deviation(p, q, u, mod)
    deviations.append((p, q, delta))

mean_dev = sum(d for _, _, d in deviations) / len(deviations)
std_dev = math.sqrt(sum((d - mean_dev) ** 2 for _, _, d in deviations) / len(deviations))

print("=== Hecke Phase Multiplicativity Analysis ===")
for (p, q, delta) in deviations:
    print(f"(p,q)=({p:3},{q:3})  Δphase = {delta:.6f} rad")
print(f"\nMean Δphase = {mean_dev:.6f} rad")
print(f"Std deviation = {std_dev:.6f} rad")

# ==============================================================
# 4. Optional: plot phase coherence
# ==============================================================

plt.figure(figsize=(8,5))
plt.hist([d for _,_,d in deviations], bins=20, color="skyblue", edgecolor="black")
plt.title("Hecke Phase Coherence Δφ Distribution")
plt.xlabel("Δφ (radians)")
plt.ylabel("Frequency")
plt.grid(True, alpha=0.3)
plt.show()