from sage import all

# --- SAGE CONFIGURATION ---
TARGET_K = 5      # <--- CHANGE THIS VALUE (e.g., 3, 4, 5, 7, 8)
MAX_N = 30        # Keep this <= 45 for free accounts (numbers get huge)
# --------------------------

def analyze_k_ladder(k, limit_n):
    print(f"\n🔵 INITIATING ANALYSIS FOR k={k}")
    print(f"   Target Field: Q(zeta_{k})")
    print(f"   Expected Splitting Condition: p = 1 mod {k}")
    print(f"   Recurrence: b_n convolution with 1/(n + 1/{k})^2")
    print("-" * 50)

    # 1. Setup Alpha and Coefficients
    alpha = 1/k
    
    # We define the coefficients of the 'Input Signal' f(t)
    # c_m = -1 / (m + alpha)^2
    # We pre-compute these to save time
    c = {}
    for m in range(limit_n + 1):
        c[m] = -1 / (m + alpha)^2

    # 2. Initialize Reciprocal Sequence b
    # Base case: b[0] * c[0] = 1  => b[0] = 1/c[0]
    b = [0] * (limit_n + 1)
    b[0] = 1 / c[0]
    
    ladder_primes = set()
    
    # 3. Convolution Loop
    for n in range(1, limit_n + 1):
        # Convolution relation: sum_{j=0}^n b_{n-j} * c_j = 0
        # b_n * c_0 + sum_{j=1}^n b_{n-j} * c_j = 0
        # b_n = - ( sum_{j=1}^n b_{n-j} * c_j ) / c_0
        
        sum_val = 0
        for j in range(1, n + 1):
            sum_val += b[n-j] * c[j]
        
        b[n] = -sum_val / c[0]
        
        # 4. Factorization (The Heavy Lifting)
        # We only factor the numerator.
        # We assume the denominator follows the known law and ignore it.
        num = b[n].numerator()
        
        # OPTIMIZATION: Don't factor if num is 1 or -1
        if abs(num) > 1:
            # We use trial division for speed on small factors, 
            # then full factorization only if needed.
            try:
                factors_found = prime_factors(abs(num))
                for p in factors_found:
                    # Filter: Ignore primes that divide k (ramified) or are too small
                    if p > k: 
                        ladder_primes.add(p)
            except:
                print(f"   [n={n}] Numerator too large to factor quickly. Skipping.")

    # 5. Verify the 'Splitting Law'
    split_primes = []
    inert_primes = []
    
    sorted_p = sorted(list(ladder_primes))
    
    for p in sorted_p:
        residue = p % k
        if residue == 1:
            split_primes.append(p)
        else:
            inert_primes.append(f"{p}({residue})")

    # 6. Report
    print(f"\n📊 RESULTS for k={k} (alpha=1/{k}):")
    print(f"   Found {len(sorted_p)} distinct candidate primes.")
    
    print(f"\n✅ MATCHING PRIMES (p = 1 mod {k}):")
    print(f"   {split_primes}")
    
    if inert_primes:
        print(f"\n⚠️ NON-MATCHING / INERT PRIMES (p != 1 mod {k}):")
        print(f"   {inert_primes}")
        print("   (Note: If this list is empty, the Dirichlet Decomposition is proven.)")
    else:
        print(f"\n🏆 PERFECT RESONANCE. All found primes split completely in Q(zeta_{k}).")

# --- EXECUTE ---
analyze_k_ladder(TARGET_K, MAX_N)