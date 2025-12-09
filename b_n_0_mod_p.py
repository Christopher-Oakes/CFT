from sage.all import *
import time

def fast_compute_bn_modp(p, N_max):
    """
    Compute b_n mod p for n=0 to N_max using power series inversion.
    This is O(N log N) per prime instead of O(N^2).
    """
    F = GF(p)
    R.<t> = PowerSeriesRing(F, default_prec=N_max+1)
    
    # Precompute C(t) = sum_{k=1}^{N_max} t^k/(2k+1)^2
    c_coeffs = [F(0)] * (N_max+1)
    for k in range(1, N_max+1):
        denom = (2*k + 1)**2
        c_coeffs[k] = F(1/denom)
    
    C = R(c_coeffs)
    D = 1 + C  # Corrected: 1 + C (was 1 - C in original)
    
    # Newton iteration for inverse: A_{n+1} = 2A_n - A_n^2 * D
    A = R(1)  # Start with A0 = 1
    current_prec = 1
    while current_prec < N_max + 1:
        next_prec = min(2 * current_prec, N_max + 1)
        A = 2*A - A*A*D
        A = A.truncate(next_prec)
        current_prec = next_prec
    
    # B(t) = (-1/4) * A(t)
    b0 = F(-1/4)
    B = b0 * A
    
    # Find zeros
    zeros = [n for n in range(N_max + 1) if B[n] == 0]
    return zeros

def find_new_repeating_primes(start_n=2001, end_n=5000, prime_bound=1000000, start_prime=37607):
    """
    Find new primes that repeat between n=start_n and n=end_n.
    """
    # Get candidate primes starting from start_prime
    all_primes = list(primes(start_prime, prime_bound))
    
    print(f"Checking {len(all_primes)} primes for repeats in n=[{start_n}, {end_n}]")
    
    repeating_primes = []
    
    for i, p in enumerate(all_primes):
        if i % 100 == 0:
            print(f"Progress: {i}/{len(all_primes)} primes checked")
        
        # Determine safe N_check
        if p > 2 * end_n + 1:
            N_check = end_n
        else:
            N_check = min(end_n, (p-1)//2 - 1)
        
        if N_check < start_n:
            continue
        
        zeros = fast_compute_bn_modp(p, N_check)
        
        # Target zeros in range
        target_zeros = [n for n in zeros if n >= start_n]
        
        # Confirm repeat
        if len(target_zeros) >= 1:
            pre_zeros = [n for n in zeros if n < start_n]
            if len(pre_zeros) >= 1 or len(target_zeros) >= 2:
                repeating_primes.append((p, sorted(pre_zeros + target_zeros)))
                print(f"New repeating prime: {p} with zeros at {sorted(pre_zeros + target_zeros)}")
    
    return repeating_primes

# RUN THE EXTENSION
if __name__ == "__main__":
    start_time = time.time()
    
    new_repeats = find_new_repeating_primes(start_n=1, end_n=1000, prime_bound=1000000, start_prime=5)
    
    end_time = time.time()
    print(f"\nFound {len(new_repeats)} new repeating primes")
    print(f"Time taken: {end_time - start_time:.2f} seconds")
    
    # Save results
    with open('new_repeating_primes_2001_5000_extended.txt', 'w') as f:
        for p, zeros in new_repeats:
            f.write(f"Prime {p}: zeros at {zeros}\n")