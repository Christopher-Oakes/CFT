import math
import json
from fractions import Fraction
from collections import defaultdict
from sympy import factorint

def compute_b_n_correct(n_max):
    """Compute coefficients b_n using the correct recurrence relation."""
    b = [Fraction(-1, 4)]  # b0
    for n in range(1, n_max + 1):
        s = Fraction(0)
        for j in range(0, n):
            denom = (2 * (n - j) + 1) ** 2
            s += b[j] / denom
        b_n = -s
        b.append(b_n)
    return b

def main():
    n_max = 2000
    print(f"Computing b_n for n=0 to {n_max} using correct recurrence...")
    b_coeffs = compute_b_n_correct(n_max)
    
    # Extract M_n for n >= 2
    print("Extracting M_n for n >= 2...")
    M_list = []  # List of tuples (n, M_n)
    for n in range(2, n_max + 1):
        numerator = b_coeffs[n].numerator
        # Check if numerator is even; if not, there's an issue
        if numerator % 2 != 0:
            print(f"Error: N_{n} is not even. Numerator: {numerator}")
            return
        M_n = abs(numerator) // 2
        M_list.append((n, M_n))
    
    print(f"Found {len(M_list)} M_n values.")
    
    # Find repeated primes using pairwise GCD
    print("Finding repeated primes via pairwise GCD...")
    prime_to_n = defaultdict(set)  # prime -> set of n values
    total_pairs = len(M_list) * (len(M_list) - 1) // 2
    processed_pairs = 0
    
    for i in range(len(M_list)):
        n_i, M_i = M_list[i]
        for j in range(i + 1, len(M_list)):
            n_j, M_j = M_list[j]
            g = math.gcd(M_i, M_j)
            if g > 1:
                factors = factorint(g)
                for prime in factors:
                    prime_to_n[prime].add(n_i)
                    prime_to_n[prime].add(n_j)
            processed_pairs += 1
            if processed_pairs % 1000 == 0:
                print(f"Processed {processed_pairs}/{total_pairs} pairs...")
    
    # Convert sets to sorted lists
    prime_to_n_sorted = {p: sorted(prime_to_n[p]) for p in prime_to_n}
    
    print(f"Found {len(prime_to_n_sorted)} primes that appear in multiple M_n values.")
    
    # Compute residues (2n+1) mod p for each prime
    print("Computing residues (2n+1) mod p...")
    prime_residues = {}
    for p, n_list in prime_to_n_sorted.items():
        residues = []
        for n in n_list:
            residue = (2 * n + 1) % p
            residues.append((n, residue))
        prime_residues[p] = residues
    
    # Save results to JSON files
    with open('repeated_prime_map_n2_2000.json', 'w') as f:
        json.dump(prime_to_n_sorted, f, indent=4)
    
    with open('repeated_prime_residues_n2_2000.json', 'w') as f:
        json.dump(prime_residues, f, indent=4)
    
    print("Analysis complete! Results saved to repeated_prime_map_n2_2000.json and repeated_prime_residues_n2_2000.json")

if __name__ == '__main__':
    main()