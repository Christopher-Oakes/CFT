import json
import math
from fractions import Fraction
from collections import defaultdict

# Compute b_n for n=0 to 300 using recurrence
def compute_b_n(max_n):
    b = [Fraction(-1, 4)]
    for n in range(1, max_n + 1):
        s = Fraction(0)
        for k in range(n):
            s += b[k] * b[n-1-k]
        denominator = 4 * (2*n + 1)**2
        b_n = -s / denominator
        b.append(b_n)
    return b

# Function to compute greatest common divisor (gcd)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Factorize a number into its prime factors
def factorize(n):
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

# Main computation
max_n = 300
print(f"Computing b_n for n=0 to {max_n}...")
b_coeffs = compute_b_n(max_n)

# Extract M_n for n>=2 (M_n = N_n / 2 if N_n is even)
print("Extracting M_n...")
M_n_dict = {}  # n: M_n
for n in range(2, max_n + 1):
    num = b_coeffs[n].numerator
    if num % 2 == 0:
        M_n_dict[n] = abs(num) // 2
    else:
        # According to the problem, for n>=2, N_n is always even, so this should not occur
        M_n_dict[n] = None

# Find repeated primes using pairwise gcd
print("Finding repeated primes via pairwise gcd...")
repeated_gcds = defaultdict(list)  # gcd value -> list of pairs (j, n)
previous_M = {}  # n: M_n for already processed n

for n in range(2, max_n + 1):
    M_n = M_n_dict.get(n)
    if M_n is None:
        continue
    for j, M_j in previous_M.items():
        g = gcd(M_n, M_j)
        if g > 1:
            repeated_gcds[g].append((j, n))
    previous_M[n] = M_n

# Factorize gcd values to get primes and map primes to n values
print("Factorizing gcd values to get primes...")
prime_to_n = defaultdict(list)  # prime p -> list of n where p divides M_n

for g, pairs in repeated_gcds.items():
    factors = factorize(g)
    primes = list(factors.keys())
    for p in primes:
        for (j, n) in pairs:
            if j not in prime_to_n[p]:
                prime_to_n[p].append(j)
            if n not in prime_to_n[p]:
                prime_to_n[p].append(n)

# Sort the n lists for each prime
for p in prime_to_n:
    prime_to_n[p].sort()

# Compute residues (2n+1) mod p for each prime and n
print("Computing residues...")
prime_to_residues = {}
for p, n_list in prime_to_n.items():
    residues = []
    for n in n_list:
        residue = (2 * n + 1) % p
        residues.append((n, residue))
    prime_to_residues[p] = residues

# Save prime_to_n as repeated_prime_map_n2_300.json
with open('repeated_prime_map_n2_300.json', 'w') as f:
    json.dump(prime_to_n, f, indent=4)

# Save prime_to_residues as repeated_prime_residues_n2_300.json
with open('repeated_prime_residues_n2_300.json', 'w') as f:
    json.dump(prime_to_residues, f, indent=4)

print("Done! Files saved: repeated_prime_map_n2_300.json and repeated_prime_residues_n2_300.json")