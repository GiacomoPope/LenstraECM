# Lenstra ECM

Implmentation of Lenstra's ECM factorisation algorithm following Pomerance's "Prime Numbers: A Computational Perspective", Algorithm 7.4.4 (Inversionless ECM).

Keen to improve this, so if you see something I should be doing differently throw some references my way.

The library relies on `gmpy2` to make modular arithmetic "fast", and a `primesieve` package for finding primes. The latter could probably be reasonably switched out with a pure python, or numpy prime sieve, but this is fast and I had it pip installed already from my work on Class Groups.

## ECM Factoring

:construction: Under construction :construction:

Maybe I could use the nice Latex support GitHub now has to give an overview of how this algorithm works... :eyes:

Calling `factor(N)` from `factor.py` first looks for small prime factors (`trial_division.py`) and then uses ECM factoring for remaining factors. The interesting code is in `lenstra_ecm.py`. Work should be done to better handle when `lenstra_ecm(N)` finds no factor. At the moment when no factor is found it tries one more time with larger bounds, then just gives up. Not so stylish...

If you wanted to use this yourself, then `import factor from factor` should just work, and the output is a dictionary of primes and exponents.

## cProfile

Running `cProfile` very clearly shows that the bulk of the computation happens in performing `xDBLADD`. Below is the output when factoring the 130 digit integer

```py
Factoring N = 2111990316278333870462901186063842291610038633220240282201156757659515273649844252195432054761834811605656662896416799950570200577
3.1923599243164062 seconds...
[(11, 1), (44756014201, 1), (152491468939, 1), (913899456083, 1), (151534339807, 1), (978212965549, 1), (39431643547271, 1), (392888992271, 1), (106325369477, 1), (240123516443, 1), (534771852913, 1), (981758150279, 1)]
```

```py
         933213 function calls (933212 primitive calls) in 3.650 seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    3.650    3.650 <string>:1(<module>)
      2/1    0.002    0.001    3.650    3.650 factor.py:5(factor)
       40    0.000    0.000    0.001    0.000 lenstra_ecm.py:30(find_curve)
       10    1.818    0.182    3.624    0.362 lenstra_ecm.py:52(lenstra_ecm)
    53278    0.094    0.000    0.094    0.000 montgomery_xz.py:1(xADD)
    53536    0.066    0.000    0.066    0.000 montgomery_xz.py:26(xDBL)
   610564    1.386    0.000    1.386    0.000 montgomery_xz.py:50(xDBLADD)
    53468    0.237    0.000    1.693    0.000 montgomery_xz.py:87(xMUL)
       40    0.000    0.000    0.000    0.000 random.py:239(_randbelow_with_getrandbits)
       40    0.000    0.000    0.000    0.000 random.py:292(randrange)
       40    0.000    0.000    0.000    0.000 random.py:366(randint)
        2    0.022    0.011    0.024    0.012 trial_division.py:4(trial_division)
     5200    0.003    0.000    0.003    0.000 {built-in method _bisect.bisect_left}
     5210    0.003    0.000    0.003    0.000 {built-in method _bisect.bisect_right}
      120    0.000    0.000    0.000    0.000 {built-in method _operator.index}
        1    0.000    0.000    3.650    3.650 {built-in method builtins.exec}
        2    0.000    0.000    0.000    0.000 {built-in method builtins.min}
       80    0.000    0.000    0.000    0.000 {built-in method builtins.pow}
        3    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.ceil}
       74    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.gcd}
       22    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.is_prime}
        3    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.isqrt}
       10    0.000    0.000    0.000    0.000 {built-in method gmpy2.gmpy2.sqrt}
    13350    0.002    0.000    0.002    0.000 {built-in method math.log}
       12    0.006    0.001    0.006    0.001 {built-in method primesieve._primesieve.primes}
    84546    0.006    0.000    0.006    0.000 {method 'append' of 'list' objects}
    53508    0.004    0.000    0.004    0.000 {method 'bit_length' of 'int' objects}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
       50    0.000    0.000    0.000    0.000 {method 'getrandbits' of '_random.Random' objects}
```

Thoughts:

- Would implementing Montgomery Arithmetic for `xDBLADD` speed it up, or slow it down due to the overhead? 11 Multiplications isn't so much...
- Zimmerman has a very different approach for Stage Two, which *also* uses affine coordinates for Stage Two. I need to learn how this works, implement and then compare. Currently Stage One is much shorter than Stage Two, and from the literature I think bounds should be balanced so they take approximately the same time?
- Would be nice to include a $P-1$ factoring if ECM fails, just incase we get lucky
- I currently don't have any check that $N \neq x^k$. Perfect powers could be tested after trial division and before ECM?
- Rather than implementing all of this in Rust/C, would there be too much FFI overhead to shift the group operations into C and then hook these back into Python?


## GMP-ECM

GMP-ECM is a highly optimised implementation of Lenstra's ECM factoring algorithm written by Paul Zimmerman and friends. It should without doubt be used instead of what I have written myself.

You can view the code in the [Inria GitLab](https://gitlab.inria.fr/zimmerma/ecm) and read a wonderful document about how it works here: [20 years of ECM, Paul Zimmermann](https://hal.inria.fr/inria-00070192v1/document). 

You can either build the binary yourself, or if you have SageMath to hand, there's bindings for it already!

```py
sage: N = 21119903162783338704629011860638422916100386332202402822011567576595152736498442521954320547618348116056566628
....: 96416799950570200577
sage: time ecm.factor(N)                                                                                                
CPU times: user 12.5 ms, sys: 49.5 ms, total: 62.1 ms
Wall time: 172 ms
[11,
 44756014201,
 106325369477,
 151534339807,
 152491468939,
 240123516443,
 392888992271,
 534771852913,
 913899456083,
 978212965549,
 981758150279,
 39431643547271]
 ```

 Then, if you're me, write a README.md for your factoring project and try not to compare the CPU time of GMP-ECM with your own implementation... Maybe, just blame python for being slow.
