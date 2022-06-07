from gmpy2 import mpz, gcd, ceil, sqrt
from math import log
from primesieve import primes
from random import randint
from bisect import bisect_left, bisect_right
from montgomery_xz import *

"""
From GMP ECM, but maybe doesn't apply
here...  
https://gitlab.inria.fr/zimmerma/ecm/

In the form:

digits : [B1, num_curves]
"""
OPTIMAL_PARAMETERS = {
    20 : [11000,     74],
    25 : [50000,     221],
    30 : [250000,    453],
    35 : [1000000,   984],
    40 : [3000000,   2541],
    45 : [11000000,  4949],
    50 : [43000000,  8266],
    55 : [110000000, 20158],
    60 : [260000000, 47173],
    65 : [850000000, 77666]
}

def find_curve(N):
    """
    Pomerance "Prime Numbers"
        Theorem 7.4.3 (ECM curve construction)

    Using a random seed, sigma
    find a Mont. curve:

    y^2 = x (x^2 + Ax + 1)
    
    With a point (maybe on the twist)

    (X : Z) = (u^3, v^3)
    """
    sigma = randint(6, N-1)
    u = (sigma**2 - 5) % N
    v = (4 * sigma)    % N
    A = (v - u)**3 * (3*u + v) * pow(4*u**3*v, -1, N) - 2
    A = A % N
    A24 = (A + 2) * pow(4, -1, N) % N
    return (N, A, A24), (u**3 % N, v**3 % N)

def lenstra_ecm(N, factor_digits=20):
    """
    Pomerance "Prime Numbers"
        Algorithm 7.4.4 (Inversionless ECM)

    B1 bounds taken from GMP-ECM
    B2 bounds taken by Pomerance (100*B1)
    D  memory bound guessed at sqrt(B2) ??
    """
    B1, num_curves = OPTIMAL_PARAMETERS[factor_digits]
    B2 = 100*B1
    # Not sure how to pick this.
    # Pomerance's guidance is D = 100
    # Other impl use sqrt of B2
    D = int(sqrt(B2))

    # Sieve for primes smaller than B2 + 2*D
    primes_under_B2  = primes(B2 + 2*D)
    # Use this to also grab primes smaller than B1
    primes_under_B1 = primes_under_B2[:bisect_right(primes_under_B2, B1)]
    
    # Precompute prime powers for stage 1
    # Used to make the power smooth integer
    # k = Π p_i**a_i; Q = [k]Q
    stage_one_prime_powers = []
    for pi in primes_under_B1:
        # pi**ai <= B1 < pi**(ai+1)
        ai = int(log(B1, pi))
        stage_one_prime_powers.append(pi**ai)

    # Precompute prime windows for stage 2
    stage_two_primes = {}
    for r in range(B1 - 1, B2, 2*D):
        r_lower_bound = r + 2
        r_upper_bound = r + 2*D
        stage_two_primes[r] = primes_under_B2[bisect_left(primes_under_B2, r_lower_bound):bisect_right(primes_under_B2, r_upper_bound)]

    for _ in range(num_curves):
        # Stage One
        # Multiply the found point Q by a
        # B1 power-smooth integer 
        E, Q = find_curve(N)        
        for pa in stage_one_prime_powers:
            Q = xMUL(Q, pa, E)
        # Lucky? if 1 < gcd(Q.Z, N) < N
        # we have a factor
        g = gcd(Q[1], N)
        if 1 < g < N: return g

        # Stage Two
        # Let's look for a factor B1 < p < B2
        S1 = xDBL(Q,  E)
        S2 = xDBL(S1, E)
        β1 = (S1[0] * S1[1]) % N
        β2 = (S2[0] * S2[1]) % N

        # Pomerance uses index one notation...
        Sd = [None, S1, S2]
        βd = [None, β1, β2]

        for d in range(3, D+1):
            S = xADD(Sd[d-1], Sd[1], Sd[d-2], E)
            Sd.append(S)
            βd.append((S[0] * S[1]) % N)

        g, B = 1, B1 - 1
        T = xMUL(Q, B - 2*D, E)
        R = xMUL(Q, B, E)

        # Dict values taken from range(B, B2, 2*D)
        for r in stage_two_primes:
            α = (R[0] * R[1] % N)
            for q in stage_two_primes[r]:
                δ = (q - r) // 2
                g = g*((R[0] - Sd[δ][0])*(R[1] + Sd[δ][1]) - α + βd[δ]) % N
            R, T = xADD(R, Sd[D], T, E), R
        g = gcd(g, N)
        if 1 < g < N: return g        

    return 0

