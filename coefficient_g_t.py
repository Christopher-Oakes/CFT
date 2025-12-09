import sympy as sp

# Define symbols
t = sp.symbols('t')
n = sp.symbols('n')

# Define f(t) expansion
N = 50  # number of terms
f_series = -4 * sp.summation(t**n / (2*n+1)**2, (n, 0, N))

# Invert series: G(t) = 1/f(t)
G_series = sp.series(1/f_series, t, 0, N+1).removeO()

# Extract coefficients
coeffs = [sp.simplify(G_series.expand().coeff(t, n)) for n in range(N+1)]

print("First coefficients of G(t) = 1/f(t):")
for n, c in enumerate(coeffs):
    print(f"b_{n} = {c}")

# Factor denominators
print("\nPrime factorizations of denominators:")
for n, c in enumerate(coeffs):
    if c.is_Rational:
        num, den = sp.fraction(c)
        print(f"b_{n}: denominator {den} -> {sp.factorint(den)}")
