import csv
import json
import math
from collections import defaultdict
from sympy import factorint

def analyze_prime_ladder(csv_filepath):
    """
    Analyze the prime-ladder phenomenon from a CSV file of coefficients.
    
    Args:
        csv_filepath (str): Path to the CSV file containing coefficients
    """
    # Read the CSV file
    data = []
    with open(csv_filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            n = int(row['n'])
            if n >= 2:  # Only consider n >= 2
                try:
                    N_n = int(row['N_n'])
                    M_n = N_n // 2  # M_n = N_n / 2
                    data.append((n, M_n))
                except (ValueError, KeyError):
                    continue
    
    print(f"Loaded {len(data)} M_n values")
    
    # Find repeated primes using pairwise GCD
    prime_to_n = defaultdict(set)
    
    # Perform pairwise GCD checks
    for i in range(len(data)):
        n_i, M_i = data[i]
        for j in range(i + 1, len(data)):
            n_j, M_j = data[j]
            g = math.gcd(M_i, M_j)
            if g > 1:
                # Factorize the GCD to get primes
                factors = factorint(g)
                for prime in factors:
                    prime_to_n[prime].add(n_i)
                    prime_to_n[prime].add(n_j)
    
    # Convert sets to sorted lists
    prime_to_n_sorted = {p: sorted(prime_to_n[p]) for p in prime_to_n}
    
    print(f"Found {len(prime_to_n_sorted)} primes that appear in multiple M_n values")
    
    # Compute residues (2n+1) mod p for each prime
    prime_residues = {}
    for p, n_list in prime_to_n_sorted.items():
        residues = []
        for n in n_list:
            residue = (2 * n + 1) % p
            residues.append((n, residue))
        prime_residues[p] = residues
    
    # Save results to JSON files
    with open('repeated_prime_map.json', 'w') as f:
        json.dump(prime_to_n_sorted, f, indent=4)
    
    with open('repeated_prime_residues.json', 'w') as f:
        json.dump(prime_residues, f, indent=4)
    
    print("Analysis complete!")
    print("Results saved to:")
    print("- repeated_prime_map.json")
    print("- repeated_prime_residues.json")
    
    # Print a summary of the results
    print("\nSummary of repeated primes:")
    for p, n_list in prime_to_n_sorted.items():
        print(f"Prime {p} appears in M_n for n = {n_list}")

if __name__ == '__main__':
    csv_path = r"C:\Users\chrys\Downloads\b_n_data_up_to_100.csv"
    analyze_prime_ladder(csv_path)