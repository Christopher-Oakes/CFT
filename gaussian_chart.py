import matplotlib.pyplot as plt
import numpy as np
from sage.all import primes, sqrt

# Your actual ladder primes under 5000
ladder_primes = [73, 103, 107, 199, 317, 421, 433, 509, 521, 547, 613, 683, 821, 823, 859, 881, 997, 1019, 1031, 1091, 1117, 1399, 1481, 1543, 1657, 1723, 1741, 1847, 1973, 2053, 2069, 2087, 2131, 2297, 2383, 2741, 2753, 2879, 2897, 3067, 3221, 3259, 3313, 3461, 3533, 3581, 3583, 3623, 3733, 3863, 3889, 3919, 4001, 4057, 4243, 4813, 4987, 5011, 5021, 5087, 5273, 5279, 5333, 5701, 5783, 5851, 5861, 6151, 6211, 6323, 6427, 6491, 6551, 6571, 6737, 6791, 7481, 7507, 7691, 7927, 8011, 8101, 8111, 8167, 8311, 8563, 8681, 8699, 9091, 9151, 9239, 9257, 9311, 9349, 9521, 9539, 9551, 9719, 9769, 9857, 9871, 9929, 10343, 10429, 10513, 10529, 10651, 10709, 10789, 10867, 11003, 11083, 11173, 11213, 11483, 11699, 11717, 11927, 12073, 12119, 12289, 12329, 12373, 12491, 12577, 13187, 13241, 13553, 13693, 13763, 13873, 14221, 14321, 14347, 14549, 14621, 14767, 14851, 15329, 15607, 15733, 16249, 16729, 16843, 17209, 17707, 18257, 19819, 20063, 20717, 20807, 21647, 22397, 22613, 23117, 23209, 24103, 24133, 24859, 25373, 25537, 25933, 26393, 27059, 27329, 27583, 28219, 29411, 30493, 32117, 33911, 35461, 35593, 36947, 37013, 67403, 91283, 124489, 2505311]


print("GENERATING GEOMETRIC PLOTS...")

# PLOT 1: Gaussian angles for p ≡ 1 mod 4 primes
def plot_gaussian_angles():
    """Plot the angles of ladder primes in complex plane"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Get p ≡ 1 mod 4 ladder primes and their Gaussian decompositions
    angles = []
    radii = []
    points = []
    
    for p in ladder_primes:
        if p % 4 == 1:
            for a in range(1, int(sqrt(p)) + 1):
                b2 = p - a*a
                if b2 >= 0:
                    b = sqrt(b2)
                    if b.is_integer():
                        b = int(b)
                        if b > a:  # Your geometric constraint
                            angle = np.arctan2(b, a)
                            angles.append(angle)
                            radii.append(np.sqrt(a*a + b*b))
                            points.append((a, b, p))
                            break
    
    # Convert to degrees for plotting
    angles_deg = [a * 180 / np.pi for a in angles]
    
    # Plot 1: Polar plot of angles
    ax1.scatter(angles_deg, radii, alpha=0.7, s=50)
    ax1.set_xlabel('Angle (degrees)')
    ax1.set_ylabel('Distance from origin')
    ax1.set_title('Ladder Primes: Gaussian Angles (p ≡ 1 mod 4)')
    ax1.grid(True)
    
    # Add sector lines
    ax1.axvline(x=45, color='red', linestyle='--', alpha=0.7, label='45° boundary')
    ax1.axvline(x=90, color='red', linestyle='--', alpha=0.7, label='90° boundary')
    ax1.legend()
    
    # Plot 2: Histogram of angles
    ax2.hist(angles_deg, bins=20, alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Angle (degrees)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Gaussian Angles')
    ax2.axvline(x=45, color='red', linestyle='--', alpha=0.7)
    ax2.axvline(x=90, color='red', linestyle='--', alpha=0.7)
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('ladder_primes_angles.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"Plotted {len(angles)} p ≡ 1 mod 4 ladder primes")
    print(f"Angle range: {min(angles_deg):.1f}° to {max(angles_deg):.1f}°")
    print(f"Mean angle: {np.mean(angles_deg):.1f}°")

# PLOT 2: Mod 60 residue classes
def plot_mod60_residues():
    """Plot the residue classes mod 60"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Get residues mod 60
    residues = [p % 60 for p in ladder_primes]
    
    # Plot 1: Bar chart of residues
    residue_counts = {}
    for r in residues:
        residue_counts[r] = residue_counts.get(r, 0) + 1
    
    ax1.bar(residue_counts.keys(), residue_counts.values(), alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Residue mod 60')
    ax1.set_ylabel('Count')
    ax1.set_title('Ladder Primes: Distribution mod 60')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Circular plot (Frobenius classes)
    angles_circle = [2 * np.pi * r / 60 for r in residues]
    radii_circle = [1] * len(residues)
    
    ax2.scatter(np.cos(angles_circle), np.sin(angles_circle), s=50, alpha=0.7)
    
    # Add circle and labels for key residues
    for r in [13, 43, 47, 19, 17, 41, 7, 23, 37, 53]:
        angle = 2 * np.pi * r / 60
        ax2.text(1.2 * np.cos(angle), 1.2 * np.sin(angle), str(r), 
                ha='center', va='center', fontsize=8)
    
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1.5, 1.5)
    ax2.set_aspect('equal')
    ax2.set_title('Frobenius Classes on mod 60 Circle')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('ladder_primes_mod60.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"Key residue classes: {sorted(set(residues))}")

# PLOT 3: Complex plane visualization with extended red line
def plot_complex_plane():
    """Plot ladder primes in complex plane"""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Find the maximum values to extend the red line appropriately
    max_coord = 0
    points_data = []
    
    for p in ladder_primes:
        if p % 4 == 1:
            for a in range(1, int(sqrt(p)) + 1):
                b2 = p - a*a
                if b2 >= 0:
                    b = sqrt(b2)
                    if b.is_integer():
                        b = int(b)
                        if b > a:
                            points_data.append((a, b, p))
                            max_coord = max(max_coord, a, b)
                            break
    
    # Plot the Gaussian integers a + bi
    for a, b, p in points_data:
        ax.scatter(a, b, s=50, alpha=0.7)
        ax.text(a, b, f' {p}', fontsize=8, alpha=0.8)
    
    # Add extended sector boundaries
    max_val = max_coord * 1.1  # Extend 10% beyond the max point
    x = np.linspace(0, max_val, 100)
    ax.plot(x, x, 'r--', alpha=0.7, linewidth=2, label='y = x (45° line)')
    
    # Add the forbidden sector shading
    x_sector = np.linspace(0, max_val, 50)
    y_sector = x_sector  # 45° line
    ax.fill_between(x_sector, y_sector, max_val, alpha=0.1, color='red', label='Forbidden sector (b ≤ a)')
    
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    ax.set_xlabel('Real part (a)')
    ax.set_ylabel('Imaginary part (b)') 
    ax.set_title('Ladder Primes in Complex Plane (p ≡ 1 mod 4)\nGaussian integers a+bi with b > a')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    
    plt.tight_layout()
    plt.savefig('ladder_primes_complex.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"Plotted {len(points_data)} Gaussian integer points")
    print(f"Extended red line to coordinate value: {max_val:.1f}")

# PLOT 4: Compare with Hecke character support (UPDATED FOR ALL PRIMES)
def plot_hecke_comparison():
    """Compare ALL ladder primes with full Hecke character support"""
    # Get the conductor 60 character
    chi = DirichletGroup(180)[19]  # Your proven character
    
    # Get maximum ladder prime to set our range
    max_ladder_prime = max(ladder_primes)
    
    # Get ALL primes up to the maximum ladder prime with χ(p) ≠ 0
    all_primes_range = list(primes(max_ladder_prime + 10000))  # Go a bit beyond for context
    hecke_support = [p for p in all_primes_range if p > 5 and chi(p) != 0]
    
    # Filter hecke support to only include primes ≤ max_ladder_prime for fair comparison
    hecke_support_upto_max = [p for p in hecke_support if p <= max_ladder_prime]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Linear scale
    x_hecke = list(range(len(hecke_support_upto_max)))
    y_hecke = hecke_support_upto_max
    
    x_ladder = list(range(len(ladder_primes)))
    y_ladder = ladder_primes
    
    ax1.plot(x_hecke, y_hecke, 'b-', alpha=0.3, label=f'Hecke support (χ(p) ≠ 0): {len(hecke_support_upto_max)} primes')
    ax1.plot(x_ladder, y_ladder, 'ro-', markersize=4, linewidth=2, label=f'Ladder primes: {len(ladder_primes)} primes')
    
    ax1.set_xlabel('Index')
    ax1.set_ylabel('Prime Value')
    ax1.set_title('Ladder Primes vs Full Hecke Character Support (Linear Scale)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Log scale to show the exponential growth difference
    ax2.plot(x_hecke, y_hecke, 'b-', alpha=0.3, label=f'Hecke support (χ(p) ≠ 0)')
    ax2.plot(x_ladder, y_ladder, 'ro-', markersize=4, linewidth=2, label=f'Ladder primes')
    ax2.set_yscale('log')
    ax2.set_xlabel('Index')
    ax2.set_ylabel('Prime Value (log scale)')
    ax2.set_title('Ladder Primes vs Full Hecke Character Support (Log Scale)\nShowing Super-Exponential Growth of Ladder')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hecke_vs_ladder_all_primes.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"Maximum ladder prime: {max_ladder_prime}")
    print(f"Hecke support up to max ladder prime: {len(hecke_support_upto_max)} primes")
    print(f"Ladder primes: {len(ladder_primes)} primes") 
    
    # Calculate growth rates
    if len(ladder_primes) > 1:
        ladder_growth = ladder_primes[-1] / ladder_primes[-2] if ladder_primes[-2] != 0 else 0
        hecke_growth = hecke_support_upto_max[-1] / hecke_support_upto_max[-2] if len(hecke_support_upto_max) > 1 and hecke_support_upto_max[-2] != 0 else 0
        
        print(f"Recent ladder growth factor: {ladder_growth:.2f}")
        print(f"Recent Hecke support growth factor: {hecke_growth:.2f}")
    
    # Convert to float for formatting
    selectivity = float(100 * len(ladder_primes) / len(hecke_support_upto_max))
    print(f"Overall selectivity: {selectivity:.2f}%")

# RUN ALL PLOTS
print("GENERATING PLOTS...")
plot_gaussian_angles()
plot_mod60_residues() 
plot_complex_plane()
plot_hecke_comparison()

print("\nALL PLOTS GENERATED AND SAVED AS PNG FILES")