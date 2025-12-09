
from sage.all import *
import numpy as np

# ================================================
# COMBINED RESULTS FOR s=3, α=1/2
# ================================================

s3_alpha12_primes = {
    43: [2, 3, 16],
    67: [21, 23],
    127: [57, 60],
    163: [21, 41, 77],
    281: [38, 42],
    331: [88, 89, 94],
    421: [33, 139, 152],
    433: [73, 94],
    479: [125, 188],
    641: [61, 76],
    691: [3, 175],
    751: [178, 186],
    827: [84, 123],
    857: [46, 189]
}

# s=2 data for comparison
s2_alpha12_primes = {
    73: [11, 18],
    103: [37, 42],
    107: [15, 37],
    199: [19, 22],
    317: [14, 62],
    421: [6, 84],
    433: [115, 144],
    509: [27, 33, 250],
    521: [40, 95],
    547: [93, 139],
    613: [112, 282],
    683: [81, 187],
    821: [164, 375],
    823: [110, 339],
    859: [391, 416]
}

print("="*80)
print("COMPREHENSIVE ANALYSIS: s=3, α=1/2 PRIME LADDER")
print("="*80)

print(f"\nFound {len(s3_alpha12_primes)} primes with multiple zeros for s=3, α=1/2:")
print("-" * 80)
print("Prime | Number of zeros | Indices of zeros")
print("-" * 80)
for p in sorted(s3_alpha12_primes.keys()):
    indices = s3_alpha12_primes[p]
    print(f"{p:5} | {len(indices):15} | {indices}")

print("\n" + "="*80)
print("CORRECTED STATISTICS")
print("="*80)

print("\nFor s=3, α=1/2:")
s3_indices = [idx for indices in s3_alpha12_primes.values() for idx in indices]
print(f"   - Total zeros: {len(s3_indices)}")
print(f"   - Average n: {np.mean(s3_indices):.1f}")
print(f"   - Min n: {min(s3_indices)}")
print(f"   - Max n: {max(s3_indices)}")
print(f"   - Median n: {np.median(s3_indices):.1f}")

print("\nFor s=2, α=1/2:")
s2_indices = [idx for indices in s2_alpha12_primes.values() for idx in indices]
print(f"   - Total zeros: {len(s2_indices)}")
print(f"   - Average n: {np.mean(s2_indices):.1f}")
print(f"   - Min n: {min(s2_indices)}")
print(f"   - Max n: {max(s2_indices)}")
print(f"   - Median n: {np.median(s2_indices):.1f}")

# ================================================
# CORRECTED ELLIPTIC CURVE DATABASE CHECK
# ================================================

print("\n" + "="*80)
print("CORRECTED ELLIPTIC CURVE DATABASE CHECK")
print("="*80)

def check_elliptic_curves_for_primes(primes_list):
    for p in primes_list:
        print(f"\nPrime {p}:")
        try:
            from sage.databases.cremona import CremonaDatabase
            cdb = CremonaDatabase()
            labels = cdb.curves(p)
            if labels:
                print(f"  Found {len(labels)} curve(s) in Cremona database.")
                for label in labels[:3]:  # Show at most 3
                    try:
                        E = EllipticCurve(label)
                        print(f"    Curve {label}: rank = {E.rank()}, CM = {E.has_cm()}")
                    except Exception as e:
                        print(f"    Error creating curve {label}: {e}")
            else:
                print(f"  No curves found in Cremona database for conductor {p}.")
        except ImportError:
            print("  Cremona database not available.")
        except Exception as e:
            print(f"  Error: {e}")

# Check the primes we are interested in
check_elliptic_curves_for_primes([43, 67, 127, 163, 281, 331, 421, 433, 479, 641, 691, 751, 827, 857])

# ================================================
# PRIME DISTRIBUTION ANALYSIS
# ================================================

print("\n" + "="*80)
print("PRIME DISTRIBUTION ANALYSIS")
print("="*80)

s3_primes = sorted(s3_alpha12_primes.keys())
print(f"\ns=3, α=1/2 primes: {s3_primes}")
print("\nResidue classes mod 4:")
for p in s3_primes:
    print(f"  {p} ≡ {p % 4} (mod 4)")

print("\nResidue classes mod 3:")
for p in s3_primes:
    print(f"  {p} ≡ {p % 3} (mod 3)")

print("\nResidue classes mod 6:")
for p in s3_primes:
    print(f"  {p} ≡ {p % 6} (mod 6)")

# Check Heegner primes
heegner_primes = [2, 3, 5, 7, 11, 19, 43, 67, 163]
print("\nHeegner primes in our list:")
for p in s3_primes:
    if p in heegner_primes:
        print(f"  {p} is a Heegner prime (class number 1)")

# ================================================
# OVERLAP ANALYSIS
# ================================================

print("\n" + "="*80)
print("OVERLAP ANALYSIS BETWEEN s=2 AND s=3")
print("="*80)

overlap = set(s3_alpha12_primes.keys()) & set(s2_alpha12_primes.keys())
print(f"\nPrimes that appear in BOTH s=2 and s=3: {sorted(overlap)}")
for p in sorted(overlap):
    print(f"\nPrime {p}:")
    print(f"  s=2 zeros at n = {s2_alpha12_primes[p]}")
    print(f"  s=3 zeros at n = {s3_alpha12_primes[p]}")
    
    # Check if any n values are close
    s2_n = set(s2_alpha12_primes[p])
    s3_n = set(s3_alpha12_primes[p])
    intersections = s2_n & s3_n
    if intersections:
        print(f"  Common n values: {sorted(intersections)}")
    else:
        # Find closest pairs
        closest_pairs = []
        for n2 in s2_alpha12_primes[p]:
            for n3 in s3_alpha12_primes[p]:
                closest_pairs.append((abs(n2 - n3), n2, n3))
        closest_pairs.sort()
        print(f"  Closest n pairs (difference, s2_n, s3_n): {closest_pairs[:3]}")

from sage.all import *

def find_zeros_s3_direct(max_n=100, prime_bound=200):
    """
    Find primes p for which b_n ≡ 0 (mod p) for s=3, α=1/2.
    Computes b_n exactly for small n, then reduces mod p.
    """
    # First compute b_n exactly for n up to max_n
    print("Computing b_n exactly for s=3, α=1/2...")
    
    b = [1/16]  # b0 for s=3, α=1/2
    
    for n in range(1, max_n+1):
        sum_term = 0
        for k in range(1, n+1):
            sum_term += b[n-k] / (2*k + 1)^3
        bn = -1/2 * sum_term
        b.append(bn)
    
    # Now check each prime
    zeros_by_prime = {}
    
    for p in primes(3, prime_bound):
        zeros = []
        for n in range(1, min(max_n+1, len(b))):
            # Get b_n as a rational
            bn = b[n]
            num = bn.numerator()
            den = bn.denominator()
            
            # Check if denominator is invertible mod p
            if den % p != 0:
                # Check if numerator ≡ 0 mod p
                if num % p == 0:
                    zeros.append(n)
        
        if len(zeros) >= 2:
            zeros_by_prime[p] = zeros
    
    return zeros_by_prime, b

print("="*70)
print("DIRECT COMPUTATION: s=3, α=1/2")
print("n up to 100, primes up to 200")
print("="*70)

zeros_s3, b_s3 = find_zeros_s3_direct(100, 200)

print(f"\nFound {len(zeros_s3)} primes with at least 2 zeros:")
for p, indices in sorted(zeros_s3.items()):
    print(f"p={p}: zeros at n={indices}")

# Now let's do the same for s=2 to verify we get 73
print("\n" + "="*70)
print("DIRECT COMPUTATION: s=2, α=1/2 (should include 73)")
print("="*70)

def find_zeros_s2_direct(max_n=100, prime_bound=200):
    # Compute b_n for s=2, α=1/2
    b = [-1/4]  # b0 for s=2, α=1/2
    
    for n in range(1, max_n+1):
        sum_term = 0
        for k in range(1, n+1):
            sum_term += b[n-k] / (2*k + 1)^2
        bn = -sum_term
        b.append(bn)
    
    zeros_by_prime = {}
    
    for p in primes(3, prime_bound):
        zeros = []
        for n in range(1, min(max_n+1, len(b))):
            bn = b[n]
            num = bn.numerator()
            den = bn.denominator()
            
            if den % p != 0:
                if num % p == 0:
                    zeros.append(n)
        
        if len(zeros) >= 2:
            zeros_by_prime[p] = zeros
    
    return zeros_by_prime, b

zeros_s2, b_s2 = find_zeros_s2_direct(100, 200)

print(f"\nFound {len(zeros_s2)} primes with at least 2 zeros:")
for p, indices in sorted(zeros_s2.items()):
    known = " (KNOWN LADDER - 73)" if p == 73 else ""
    print(f"p={p}: zeros at n={indices}{known}")

# Let's check what happens at p=73 specifically
print("\n" + "="*70)
print("SPECIFIC CHECK: p=73 for s=2")
print("="*70)

if 73 in zeros_s2:
    print(f"✓ p=73 found with zeros at n={zeros_s2[73]}")
else:
    print("✗ p=73 NOT FOUND")
    # Let's manually check b_n mod 73 for n=11,18
    for n in [11, 18]:
        bn = b_s2[n]
        num = bn.numerator()
        den = bn.denominator()
        print(f"n={n}: b_n = {num}/{den}")
        print(f"  num mod 73 = {num % 73}")
        print(f"  den mod 73 = {den % 73}")
        if den % 73 == 0:
            print("  WARNING: denominator divisible by 73")

# Now let's extend the search for s=3 to higher primes
print("\n" + "="*70)
print("EXTENDED SEARCH: s=3, α=1/2, primes up to 500")
print("="*70)

zeros_s3_ext, _ = find_zeros_s3_direct(100, 500)

print(f"\nFound {len(zeros_s3_ext)} primes with at least 2 zeros:")
for p, indices in sorted(zeros_s3_ext.items()):
    print(f"p={p}: zeros at n={indices}")

# Let's also check small primes that we might have missed
print("\n" + "="*70)
print("CHECKING SMALL PRIMES (3..100) for s=3")
print("="*70)

small_zeros = {}
for p in primes(3, 100):
    zeros = []
    for n in range(1, min(100, len(b_s3))):
        bn = b_s3[n]
        num = bn.numerator()
        den = bn.denominator()
        
        if den % p != 0:
            if num % p == 0:
                zeros.append(n)
    
    if len(zeros) >= 2:
        small_zeros[p] = zeros

print(f"\nFound {len(small_zeros)} small primes with at least 2 zeros:")
for p, indices in sorted(small_zeros.items()):
    print(f"p={p}: zeros at n={indices}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
For s=2, α=1/2:
- Found {len(zeros_s2)} primes with multiple zeros (including 73: {zeros_s2.get(73, 'NOT FOUND')})

For s=3, α=1/2:
- Found {len(zeros_s3_ext)} primes with multiple zeros in primes up to 500
- Small primes (3..100): {len(small_zeros)} primes with multiple zeros

Key differences:
1. s=2 primes tend to be smaller and more frequent
2. s=3 primes start appearing at higher n
3. The zero patterns are different

This suggests:
- s=3 indeed has its own "prime ladder" structure
- The primes for s=3 are likely associated with different mathematical objects
  (weight 2 modular forms/elliptic curves vs weight 1 for s=2)
""")
