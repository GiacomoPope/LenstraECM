from trial_division import trial_division
from lenstra_ecm import lenstra_ecm
from perfect_power import perfect_power
from gmpy2 import is_prime, is_power


def factor(n, perfect_power_exp=1):

    def update_factors(N, q, factors):
        while N % q == 0:
            if q in factors:
                factors[int(q)] += perfect_power_exp
            else:
                factors[int(q)] = perfect_power_exp
            N //= q
        return N, factors

    factors, N = trial_division(n)
    while True:
        # All done
        if N == 1:
            break
        # Remaining integer is prime; Done!
        elif is_prime(N):
            factors[int(N)] = perfect_power_exp
            break
        # Check if N can be written as (n)^e
        elif is_power(N):
            N, e = perfect_power(N)
            if is_prime(N):
                if N in factors:
                    factors[int(N)] += e*perfect_power_exp
                else:
                    factors[int(N)] = e*perfect_power_exp
                break
            else:
                perfect_power_exp = e*perfect_power_exp

        # Start finding factors with ECM
        q = lenstra_ecm(N)
        # We found no factor, try again... 
        if q == 0:
            print(f"No factor found, increasing bounds")
            factor_digits = int(len(str(N)) / 10) * 5
            factor_digits = max(factor_digits, 25)
            q = lenstra_ecm(N, factor_digits=factor_digits)
            if q == 0:
                print(f"No factor found, giving up...")
                print(f"Remaining composite: {N=}")
                return factors
        # ECM can return a non-prime factor
        # If this is the case, factor again!
        if not is_prime(q):
            sub_factors = factor(q, perfect_power_exp=perfect_power_exp)
            for sub_factor in sub_factors:
                N, factors = update_factors(N, sub_factor, factors)
        else:
            N, factors = update_factors(N, q, factors)
    return factors

if __name__ == '__main__':
    from gmpy2 import next_prime
    from random import seed, randint
    from factor_util import test_factors
    seed(0)

    # Generate some high power containing 
    # number to make sure new work works!

    # Factor me!
    n = randint(1, 2**30)
    n = pow(n, randint(1,5))
    e1, e2 = randint(1,10), randint(1,10)
    k = randint(1,10)
    for _ in range(5):
        p = next_prime(randint(1, 2**40))
        n *= p**e1
    for _ in range(3):
        p = next_prime(randint(1, 2**40))
        n *= p**e2
    n = pow(n,k)
    n *= next_prime(randint(1, 2**40))

    if 0:
        """
        Check the factoring
        """
        from time import time
        print(f"Factoring N = {n}")
        t = time()
        factors = factor(n)
        print(f"{time() - t} seconds...")
        assert test_factors(n, factors)
        print(sorted(factors.items()))
    else:
        """
        Bench the factoring
        """
        import cProfile
        gvars = {'n' : n}
        lvars = {'factor': factor}
        cProfile.runctx('factor(n)', globals=gvars, locals=lvars)
