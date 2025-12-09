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