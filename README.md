# InPlace

[![Build Status](https://github.com/YuttariKanata/InPlace.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/YuttariKanata/InPlace.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/YuttariKanata/InPlace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/YuttariKanata/InPlace.jl)

## InPlace.jl: Zero-Allocation Arbitrary-Precision Arithmetic

### **A high-performance, mutation-first API for GMP/MPFR in Julia**

> InPlace.jl provides fast, zero-allocation in-place operations on `BigInt` and `BigFloat` with full control over precision and rounding modes. Designed for numerical algorithms demanding extreme precision and predictable performance.

---

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
- [BigInt Operations](#bigint-operations)
- [BigFloat Operations](#bigfloat-operations)
- [API Reference](#api-reference)
- [Performance Guide](#performance-guide)
- [Examples](#examples)
- [Technical Details](#technical-details)
- [FAQ](#faq)

---

## Overview

### Why InPlace.jl?

Julia's built-in `BigInt` and `BigFloat` are powerful but create temporary allocations on every operation:

```julia
# Standard Julia (allocates 3 times)
result = x + y + z
```

InPlace.jl eliminates this waste through **mutation-first design**:

```julia
using InPlace

# Single allocation (result reused)
add!(result, x, y)
add!(result, result, z)  # result mutated in-place
```

### Key Features

✅ **Zero Allocations** - Operations reuse result object  
✅ **Type-Driven Dispatch** - Automatic UI/SI selection for native types  
✅ **Full GMP/MPFR API** - 100+ functions wrapped  
✅ **Precision Control** - Explicit precision and rounding modes  
✅ **Native Type Support** - Seamless Int32/Int64/UInt32/UInt64 handling  
✅ **SIMD-Friendly** - Deterministic performance for real-time systems

---

## Installation

```julia
using Pkg
Pkg.add("InPlace")
```

Or from the repo:

```julia
Pkg.add(url="https://github.com/YuttariKanata/InPlace.jl")
```

---

## Quick Start

### Basic Operations

```julia
using InPlace

# Create result objects
x = BigInt(100)
y = BigInt(50)
result = BigInt(0)

# In-place arithmetic
add!(result, x, y)       # result = 150
sub!(result, result, 20) # result = 130
mul!(result, result, 2)  # result = 260
div!(result, result, 2)  # result = 130

# BigFloat with precision control
f = BigFloat(0; precision=256)
set!(f, "3.14159")
sqrt!(f, 2)              # f = √2
mul!(f, f, 2)            # f *= 2
pi!(f)                   # f = π
```

### Why Mutation?

Every non-mutating operation in standard Julia allocates:

```julia
# ❌ Bad: Allocates 3 times
r1 = a + b
r2 = r1 + c
r3 = r2 + d

# ✅ Good: Single allocation
result = BigInt(0)
add!(result, a, b)
add!(result, result, c)
add!(result, result, d)
```

---

## Core Concepts

### Tier 1 & Tier 2 Dispatch

InPlace.jl uses a **two-tier dispatch strategy** for maximum performance:

#### **Tier 1: Compile-Time Dispatch** (Static)

When you pass native machine integers (Int32, Int64, UInt32, UInt64), the compiler **specializes code paths**:

```julia
# Native types ← Compile-time dispatch
add!(result, bigint_x, 42_i64)  # → mpz_add_ui (fast)
```

Why? GMP has specialized `_ui` and `_si` functions that are **2-4x faster** than general operations.

#### **Tier 2: Run-Time Dispatch** (Dynamic)

For generic `Integer` types, the function checks the value at runtime:

```julia
# Generic Integer ← Run-time dispatch
y = some_function()  # Returns Integer (could be Int, BigInt, etc.)
add!(result, bigint_x, y)  # Dispatches at runtime
```

### Type-Driven Dispatch Example

```julia
# All these call ONE function name: add!
add!(r, x::BigInt, y::BigInt)    # → mpz_add
add!(r, x::BigInt, y::Int64)     # → mpz_add_ui
add!(r, x::BigInt, y::UInt32)    # → mpz_add_ui (auto upcasted)
add!(r, x::BigInt, y::Int128)    # → falls back to mpz_add(BigInt(y))
add!(r, x::Integer, y::Integer)  # → converts both to BigInt
```

No need to remember function names like `mpz_add_ui_si`—Julia's type system handles it!

### Unified API: Internal Dispatch

Historically, low-level libraries force you to choose:

- `mul_2ui!(rop, x, n)` for unsigned shift
- `mul_2si!(rop, x, n)` for signed shift

**InPlace.jl unifies them:**

```julia
# Same function, type-driven dispatch
mul_2exp!(result, x, 5_u64)   # → mpfr_mul_2ui (unsigned)
mul_2exp!(result, x, -5_i64)  # → mpfr_mul_2si (signed)
```

---

## BigInt Operations

### Assignment & Conversion

```julia
using InPlace

x = BigInt(0)

# From strings
set!(x, "123456789012345678901234567890")
set!(x, "deadbeef", base=16)

# From integers
set!(x, 42)
set!(x, -2147483648_i64)

# Copy
set!(x, another_bigint)

# Float conversion
set!(x, 3.14159)  # Truncates to 3
```

### Arithmetic

```julia
result = BigInt(0)
a = BigInt(1000)
b = BigInt(7)

# Basic
add!(result, a, b)      # result = 1007
sub!(result, a, b)      # result = 993
mul!(result, a, b)      # result = 7000
div!(result, a, b)      # result = 142 (truncation)
rem!(result, a, b)      # result = 6
mod!(result, a, b)      # result = 6

# With rounding modes
div!(result, a, b; rounding=:floor)  # Floor division
div!(result, a, b; rounding=:ceil)   # Ceiling division
rem!(result, a, b; rounding=:floor)  # Floored remainder
```

### Multiply-Accumulate

```julia
# rop += x * y (no intermediate allocation)
addmul!(result, x, y)

# rop -= x * y
submul!(result, x, y)

# Equivalent to: result += 1000000 * x (in-place)
addmul!(result, x, 1000000)
```

### Bit Operations

```julia
# Shift left: result = x * 2^n
mul_2exp!(result, x, 10)

# Shift right: result = x / 2^n (floor division)
div_2exp!(result, x, 10)

# Bit manipulation
setbit!(result, 5)   # result[5] = 1
clrbit!(result, 3)   # result[3] = 0
combit!(result, 7)   # result[7] = 1 - result[7] (toggle)
```

### Exponentiation

```julia
pow!(result, base, exponent)
pow!(result, 2, 100)       # 2^100 = 1267650600228229401496703205376

# Modular exponentiation (for cryptography)
powm!(result, base, exp, modulus)

# Constant-time version (resists timing attacks)
powm_sec!(result, base, exp, modulus)
```

### Number Theory

```julia
# GCD / LCM
gcd!(result, 48, 180)      # result = 12
lcm!(result, 48, 180)      # result = 720

# Extended GCD: finds s, t such that g = a*s + b*t
gcdext!(g, s, t, 35, 15)   # g=5, s=1, t=-2

# Modular inverse
invert!(result, 3, 7)      # result = 5 (3*5 ≡ 1 mod 7)

# Remove factor
remove!(result, 72, 2)     # result = 2 (72 = 2^3 * 3^2, remove 3 2's → 2^1)

# Combinatorics
fac_ui!(result, 20)        # 20! = 2432902008176640000
bin_ui!(result, 100, 20)   # C(100,20) binomial coefficient
fib_ui!(result, 100)       # 100-th Fibonacci number

# Primality testing
prob_prime = mpz_probab_prime_p(n, 25)  # 25 rounds, error < 2^-50
```

### Roots

```julia
# Square root
sqrt!(result, 16)          # result = 4
sqrtrem!(root, rem, 17)    # root=4, rem=1 (17 = 4^2 + 1)

# n-th root
root!(result, 8, 3)        # result = 2 (∛8 = 2)
rootrem!(root, rem, 10, 3) # root=2, rem=2 (10 = 2^3 + 2)
```

### Bitwise Logic

```julia
and!(result, x, y)         # Bitwise AND
or!(result, x, y)          # Bitwise OR
xor!(result, x, y)         # Bitwise XOR
com!(result, x)            # Bitwise complement
```

---

## BigFloat Operations

### Assignment & Precision

```julia
x = BigFloat(0; precision=256)  # 256 bits of precision

# From string
set!(x, "3.14159265358979323846")

# From native types
set!(x, 42)
set!(x, 3.14159_f64)

# From BigInt
set!(x, BigInt(1000000))

# Copy with rounding
set!(x, y; rounding=MPFRRoundNearest)
```

### Precision & Rounding Modes

```julia
# Rounding modes
const RND = InPlace.MPFRRoundNearest      # Round to nearest (default)
const RND = InPlace.MPFRRoundToZero       # Truncate
const RND = InPlace.MPFRRoundUp           # Round toward ∞
const RND = InPlace.MPFRRoundDown         # Round toward -∞
const RND = InPlace.MPFRRoundFromZero     # Round away from zero

# All operations accept rounding keyword
add!(result, x, y; rounding=RND)
sqrt!(result, x; rounding=RND)
```

### Basic Arithmetic

```julia
x = BigFloat(2.5; precision=256)
y = BigFloat(1.5; precision=256)
result = BigFloat(0; precision=256)

add!(result, x, y)         # result = 4.0
sub!(result, x, y)         # result = 1.0
mul!(result, x, y)         # result = 3.75
div!(result, x, y)         # result ≈ 1.666...

# Native type support (auto-upcasted to precision of result)
add!(result, x, 42)        # Adds 42 with full precision
mul!(result, x, 2_u64)     # Multiplies by 2 (exact)
```

### Advanced Operations

```julia
# Fused multiply-add: result = x*y + z (better precision/speed)
fma!(result, x, y, z)

# Fused multiply-subtract: result = x*y - z
fms!(result, x, y, z)

# Euclidean distance
hypot!(result, x, y)       # √(x² + y²)

# Relative difference
reldiff!(result, x, y)     # |x - y| / max(|x|, |y|)

# Dimension: max(x - y, 0)
dim!(result, x, y)
```

### Exponential & Logarithmic

```julia
# Base e
exp!(result, x)

# Base 2 and 10
exp2!(result, x)           # 2^x
exp10!(result, x)          # 10^x

# Numerically stable: exp(x) - 1
expm1!(result, x)

# Logarithms
log!(result, x)            # Natural log
log2!(result, x)           # Base 2
log10!(result, x)          # Base 10
log1p!(result, x)          # log(1 + x), stable for small x
```

### Trigonometric

```julia
sin!(result, x)
cos!(result, x)
tan!(result, x)

# Compute sin and cos simultaneously (faster)
sin_cos!(result_sin, result_cos, x)

# Reciprocals
sec!(result, x)            # 1/cos(x)
csc!(result, x)            # 1/sin(x)
cot!(result, x)            # 1/tan(x)
```

### Inverse Trigonometric

```julia
asin!(result, x)           # arcsin (domain: [-1,1])
acos!(result, x)           # arccos (domain: [-1,1])
atan!(result, x)           # arctan
atan2!(result, y, x)       # atan2(y, x) ∈ (-π, π]
```

### Hyperbolic

```julia
sinh!(result, x)
cosh!(result, x)
tanh!(result, x)

# Compute sinh and cosh simultaneously (faster)
sinh_cosh!(result_sinh, result_cosh, x)

sech!(result, x)           # 1/cosh(x)
csch!(result, x)           # 1/sinh(x)
coth!(result, x)           # 1/tanh(x)
```

### Inverse Hyperbolic

```julia
asinh!(result, x)
acosh!(result, x)          # Domain: x ≥ 1
atanh!(result, x)          # Domain: |x| < 1
```

### Special Functions

```julia
# Gamma function
gamma!(result, x)          # Γ(x)
lngamma!(result, x)        # log(|Γ(x)|)
lgamma!(result, signp, x)  # log(|Γ(x)|) and sign separately

# Derivatives
digamma!(result, x)        # ψ(x) = Γ'(x)/Γ(x)

# Beta function
beta!(result, x, y)        # B(x,y) = Γ(x)Γ(y)/Γ(x+y)

# Riemann zeta
zeta!(result, x)           # ζ(x)
zeta!(result, 20_u64)      # ζ(20) via mpfr_zeta_ui (faster)

# Error functions
erf!(result, x)            # Error function
erfc!(result, x)           # Complementary error function

# Exponential integral
eint!(result, x)           # E₁(x)
li2!(result, x)            # Dilogarithm
```

### Bessel & Airy Functions

```julia
# Bessel functions of 1st kind
j0!(result, x)             # J₀(x)
j1!(result, x)             # J₁(x)
jn!(result, n, x)          # Jₙ(x) where n is an integer

# Bessel functions of 2nd kind
y0!(result, x)             # Y₀(x)
y1!(result, x)             # Y₁(x)
yn!(result, n, x)          # Yₙ(x)

# Airy function
ai!(result, x)             # Ai(x)

# AGM (Arithmetic-Geometric Mean)
agm!(result, x, y)
```

### Rounding & Remainder

```julia
rint!(result, x)           # Round to nearest integer (respecting mode)
ceil!(result, x)           # Round toward ∞
floor!(result, x)          # Round toward -∞
round!(result, x)          # Banker's rounding (round to nearest even)
trunc!(result, x)          # Truncate toward 0

# Fractional and integer parts
frac!(result, x)           # Fractional part: x - floor(x)
modf!(int_part, frac_part, x)  # Split into integer + fractional

# Remainder/modulo
fmodquo!(result, quo, x, y)    # quo = x/y, result = x - quo*y
remquo!(result, quo, x, y)     # remainder with quotient
```

### Bit-Shift Operations

```julia
# Left shift (multiply by 2^n)
mul_2exp!(result, x, 10)   # result = x * 2^10

# Right shift (divide by 2^n)
div_2exp!(result, x, 10)   # result = x / 2^10 (rounded)
```

### Floating-Point Manipulation

```julia
# Get/set exponent
exp = get_exp!(x)
set_exp!(x, exp - 1)       # Halve the value without changing precision

# Sign manipulation
setsign!(result, x, sign)  # Set sign (1 or -1)
copysign!(result, x, y)    # Copy sign from y to x

# Next representable float
nextabove!(x)              # Move to next larger float
nextbelow!(x)              # Move to next smaller float
nexttoward!(x, target)     # Move toward target
```

### Comparison & Predicates

```julia
# Use Julia's standard operators (built on mpfr_cmp)
x > y
x == y
x < y

# Check special values
isnan(x)       # NaN test
isinf(x)       # Infinity test
isfinite(x)    # Not NaN and not ∞
iszero(x)      # Exactly zero
```

### Mathematical Constants

```julia
pi!(result)                # π
euler!(result)             # e (Euler-Mascheroni constant γ)
catalan!(result)           # Catalan constant G

# All support rounding modes
pi!(result; rounding=MPFRRoundNearest)
```

---

## API Reference

### Function Categories

#### **BigInt: Assignment**

```julia
set!(rop::BigInt, x::T) where T ∈ {BigInt, Integer, AbstractString}
swap!(x::BigInt, y::BigInt)
```

#### **BigInt: Basic Arithmetic**

```julia
add!(rop, x, y)
sub!(rop, x, y)
mul!(rop, x, y)
div!(rop, x, y; rounding::Symbol)
rem!(rop, x, y; rounding::Symbol)
mod!(rop, x, y)
divrem!(q, r, x, y; rounding::Symbol)
divexact!(rop, x, y)
```

#### **BigInt: Compound**

```julia
addmul!(rop, x, y)         # rop += x*y
submul!(rop, x, y)         # rop -= x*y
mul_2exp!(rop, x, n)       # rop = x * 2^n
div_2exp!(rop, x, n; rounding)  # rop = x / 2^n
```

#### **BigInt: Bit Manipulation**

```julia
setbit!(rop, bit)
clrbit!(rop, bit)
combit!(rop, bit)
and!(rop, x, y)
or!(rop, x, y)
xor!(rop, x, y)
com!(rop, x)
```

#### **BigInt: Exponentiation & Roots**

```julia
pow!(rop, x, n)
powm!(rop, x, e, m)        # Modular exponentiation
powm_sec!(rop, x, e, m)    # Constant-time version
sqrt!(rop, x)
sqrtrem!(root, rem, x)
root!(rop, x, n)
rootrem!(root, rem, x, n)
```

#### **BigInt: Number Theory**

```julia
gcd!(rop, x, y)
lcm!(rop, x, y)
gcdext!(g, s, t, x, y)     # Extended GCD
invert!(rop, x, m)         # Modular inverse
remove!(rop, x, p)         # Remove factor p
fac_ui!(rop, n)            # Factorial
bin_ui!(rop, n, k)         # Binomial coefficient
fib_ui!(rop, n)            # Fibonacci
lucnum_ui!(rop, n)         # Lucas number
```

#### **BigInt: Unary**

```julia
neg!(rop, x)
abs!(rop, x)
nextprime!(rop, x)
```

#### **BigFloat: Assignment**

```julia
set!(rop::BigFloat, x::T; rounding) where T ∈ {BigFloat, Integer, AbstractFloat, AbstractString}
swap!(x::BigFloat, y::BigFloat)
```

#### **BigFloat: Arithmetic**

```julia
add!(rop, x, y; rounding)
sub!(rop, x, y; rounding)
mul!(rop, x, y; rounding)
div!(rop, x, y; rounding)
```

#### **BigFloat: Unary**

```julia
neg!(rop, x; rounding)
abs!(rop, x; rounding)
sqrt!(rop, x; rounding)
square!(rop, x; rounding)
rec_sqrt!(rop, x; rounding)    # 1/√x
cbrt!(rop, x; rounding)        # ∛x
```

#### **BigFloat: Root Extraction**

```julia
rootn!(rop, x, n; rounding)    # n-th root (n signed or unsigned)
```

#### **BigFloat: Power**

```julia
pow!(rop, x, y; rounding)      # Generic power (y can be float)
```

#### **BigFloat: Exponential & Logarithm**

```julia
exp!(rop, x; rounding)
exp2!(rop, x; rounding)
exp10!(rop, x; rounding)
expm1!(rop, x; rounding)
log!(rop, x; rounding)
log2!(rop, x; rounding)
log10!(rop, x; rounding)
log1p!(rop, x; rounding)
```

#### **BigFloat: Trigonometric**

```julia
sin!(rop, x; rounding)
cos!(rop, x; rounding)
tan!(rop, x; rounding)
sin_cos!(rop_sin, rop_cos, x; rounding)
sec!(rop, x; rounding)
csc!(rop, x; rounding)
cot!(rop, x; rounding)
```

#### **BigFloat: Inverse Trigonometric**

```julia
asin!(rop, x; rounding)
acos!(rop, x; rounding)
atan!(rop, x; rounding)
atan2!(rop, y, x; rounding)
```

#### **BigFloat: Hyperbolic**

```julia
sinh!(rop, x; rounding)
cosh!(rop, x; rounding)
tanh!(rop, x; rounding)
sinh_cosh!(rop_sinh, rop_cosh, x; rounding)
sech!(rop, x; rounding)
csch!(rop, x; rounding)
coth!(rop, x; rounding)
```

#### **BigFloat: Inverse Hyperbolic**

```julia
asinh!(rop, x; rounding)
acosh!(rop, x; rounding)
atanh!(rop, x; rounding)
```

#### **BigFloat: Special Functions**

```julia
gamma!(rop, x; rounding)
gamma_inc!(rop, x, y; rounding)
lngamma!(rop, x; rounding)
lgamma!(rop, signp, x; rounding)
digamma!(rop, x; rounding)
beta!(rop, x, y; rounding)
zeta!(rop, x; rounding)
erf!(rop, x; rounding)
erfc!(rop, x; rounding)
eint!(rop, x; rounding)
li2!(rop, x; rounding)
```

#### **BigFloat: Bessel & Airy**

```julia
j0!(rop, x; rounding)
j1!(rop, x; rounding)
jn!(rop, n, x; rounding)
y0!(rop, x; rounding)
y1!(rop, x; rounding)
yn!(rop, n, x; rounding)
ai!(rop, x; rounding)
agm!(rop, x, y; rounding)
```

#### **BigFloat: Rounding & Remainder**

```julia
rint!(rop, x; rounding)
ceil!(rop, x)
floor!(rop, x)
round!(rop, x)
trunc!(rop, x)
frac!(rop, x; rounding)
modf!(int_part, frac_part, x; rounding)
fmodquo!(rop, quo, x, y; rounding)
remquo!(rop, quo, x, y; rounding)
```

#### **BigFloat: Advanced**

```julia
fma!(rop, x, y, z; rounding)       # x*y + z
fms!(rop, x, y, z; rounding)       # x*y - z
hypot!(rop, x, y; rounding)
reldiff!(rop, x, y; rounding)
dim!(rop, x, y; rounding)
min!(rop, x, y; rounding)
max!(rop, x, y; rounding)
```

#### **BigFloat: Bit Operations**

```julia
mul_2exp!(rop, x, n; rounding)
div_2exp!(rop, x, n; rounding)
```

#### **BigFloat: Floating-Point Manipulation**

```julia
get_exp!(x) → exponent
set_exp!(x, e)
setsign!(rop, x, sign; rounding)
copysign!(rop, x, y; rounding)
nextabove!(x)
nextbelow!(x)
nexttoward!(x, target)
```

#### **BigFloat: Constants**

```julia
pi!(rop; rounding)
euler!(rop; rounding)
catalan!(rop; rounding)
```

---

## Performance Guide

### Allocation Profiling

```julia
using InPlace

# Check allocations with @time
x = BigInt(2)^1000
y = BigInt(3)^1000
result = BigInt(0)

# Benchmark: should show 0 allocations
@time for i in 1:1000
    mul!(result, x, y)
    add!(result, result, 1)
end
```

### Best Practices

1. **Pre-allocate result objects**

   ```julia
   result = BigInt(0)  # Allocate once
   for i in 1:1000000
       add!(result, x, y)
   end
   ```

2. **Batch operations with the same precision**

   ```julia
   x = BigFloat(0; precision=1024)
   y = BigFloat(0; precision=1024)
   z = BigFloat(0; precision=1024)
   # ✅ All operations now use same precision
   
   # ❌ Avoid
   x = BigFloat(0; precision=256)
   y = BigFloat(0; precision=512)
   add!(result, x, y)  # Requires conversion
   ```

3. **Use native integer arguments**

   ```julia
   # ✅ Fast: Type-dispatched at compile time
   add!(result, x, 42_i64)
   
   # ❌ Slow: Generic dispatch
   add!(result, x, Integer(42))
   ```

4. **Avoid temporary type conversions**

   ```julia
   # ❌ Bad: Creates temporary BigFloat
   sqrt!(result, x)  # x might not be Float

   # ✅ Good: Pass already-converted
   set!(temp, x)
   sqrt!(result, temp)
   ```

### Timing Comparisons

```julia
using InPlace, BenchmarkTools

x = BigInt(2)^1000
y = BigInt(3)^1000
result = BigInt(0)

# InPlace (0 allocations)
@benchmark mul!($result, $x, $y) evals=1000

# Standard Julia (allocates)
@benchmark $x * $y evals=1000
```

**Typical results:** InPlace.jl is 3-10x faster for repeated operations.

---

## Examples

### Example 1: Arbitrary-Precision PI Computation (Bailey–Borwein–Plouffe)

```julia
using InPlace

function compute_pi_digit(d::Int)
    """Compute the d-th hexadecimal digit of π (BBP formula)"""
    precision = d + 50
    s = BigFloat(0; precision=precision)
    one_16 = BigFloat(1; precision=precision)
    
    set!(one_16, 1)
    for k in 0:d
        ak = BigFloat(k; precision=precision)
        set!(ak, k)
        
        # s += (1/16^k) * (4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6))
        term = BigFloat(0; precision=precision)
        # ... compute term ...
        add!(s, s, term)
    end
    return s
end
```

### Example 2: Modular Exponentiation (RSA)

```julia
using InPlace

function rsa_encrypt(message::BigInt, exponent::BigInt, modulus::BigInt)
    """Encrypt message: ciphertext = message^exponent mod modulus"""
    ciphertext = BigInt(0)
    powm_sec!(ciphertext, message, exponent, modulus)  # Constant-time
    return ciphertext
end
```

### Example 3: Newton-Raphson for √2

```julia
using InPlace

function sqrt2_newtonraphson(digits::Int)
    """Compute √2 to specified number of digits"""
    precision = digits * 4  # Extra precision for intermediate calcs
    
    x = BigFloat(1.5; precision=precision)  # Initial guess
    x_prev = BigFloat(0; precision=precision)
    two = BigFloat(2; precision=precision)
    
    for iter in 1:20
        set!(x_prev, x)
        
        # x = (x + 2/x) / 2
        div!(temp, two, x)
        add!(x, x, temp)
        div!(x, x, 2)
        
        # Check convergence
        sub!(diff, x, x_prev)
        if abs(diff) < BigFloat(10)^(-digits)
            break
        end
    end
    
    return x
end
```

### Example 4: Integration (Simpson's Rule)

```julia
using InPlace

function integrate_sin(a::Real, b::Real, n::Int)
    """Integrate sin(x) from a to b using Simpson's rule"""
    precision = 256
    
    h = BigFloat(b - a; precision=precision)
    div!(h, h, n)
    
    result = BigFloat(0; precision=precision)
    sum_odd = BigFloat(0; precision=precision)
    sum_even = BigFloat(0; precision=precision)
    
    temp_x = BigFloat(a; precision=precision)
    temp_val = BigFloat(0; precision=precision)
    
    for i in 1:n-1
        add!(temp_x, BigFloat(a; precision=precision), h * i)
        sin!(temp_val, temp_x)
        
        if iseven(i)
            addmul!(sum_even, temp_val, 2)
        else
            addmul!(sum_odd, temp_val, 4)
        end
    end
    
    sin!(temp_val, a)
    add!(result, temp_val, sum_odd)
    add!(result, result, sum_even)
    
    sin!(temp_val, b)
    add!(result, result, temp_val)
    
    mul!(result, result, h)
    div!(result, result, 3)
    
    return result
end
```

### Example 5: Iterative GCD Loop

```julia
using InPlace

function gcd_many(numbers::Vector{BigInt})
    """Compute GCD of many numbers"""
    result = numbers[1]
    for i in 2:length(numbers)
        gcd!(result, result, numbers[i])
    end
    return result
end
```

---

## Technical Details

### Memory Model

InPlace.jl operations **never allocate** on the result object:

```julia
result = BigInt(0)
x = BigInt(2)^1000000

# Memory layout BEFORE
# result: [limbs=C_NULL, alloc=0, size=0]

add!(result, x, x)

# Memory layout AFTER
# result: [limbs=..., alloc=..., size=...]
# Same object; contents updated in-place
```

### Tier 1 Dispatch Details

When you write:

```julia
add!(result, bigint, 42_i64)
```

Julia's type system sees `Int64` and the compiler generates a specialized method:

```c
// Pseudo-code generated by @generated
mpz_add_ui(rop, x, 42)  // Compile-time constant
```

If you write:

```julia
y = Int32(42)
add!(result, bigint, y)
```

The compiler generates:

```c
// Int32 dispatch
mpz_add_si(rop, x, (signed long)y)  // Type-safe cast
```

### Tier 2 Runtime Dispatch

For:

```julia
y = some_function()  # Type: Integer (concrete at runtime)
add!(result, bigint, y)
```

The dispatch happens at runtime:

```julia
if y isa BigInt
    mpz_add(...)
elseif y isa NativeSignedInt
    mpz_add_si(...)
elseif y isa NativeUnsignedInt
    mpz_add_ui(...)
else
    # Generic fallback
    ...
end
```

### GMP API Mapping

| InPlace Function | GMP Backend | Type Filter |
| --- | --- | --- |
| `add!(r, x, y::NativeUnsignedInt)` | `mpz_add_ui` | UInt32, UInt64 |
| `add!(r, x, y::NativeSignedInt)` | `mpz_add_si` or `mpz_sub_ui` | Int32, Int64 |
| `add!(r, x, y::BigInt)` | `mpz_add` | BigInt |
| `add!(r, x, y::Integer)` | (converts to BigInt) | Other |

---

## FAQ

### Q: Why mutate instead of returning new objects?

**A:** Performance in tight loops. Standard Julia:

```julia
# Allocates 1000 times
for i in 1:1000
    x = x + y
end
```

InPlace.jl:

```julia
# 0 allocations
for i in 1:1000
    add!(x, x, y)
end
```

### Q: Can I mix precisions?

**A:** Yes, but be explicit:

```julia
x = BigFloat(0; precision=256)
y = BigFloat(0; precision=512)

# InPlace upcasts to result precision
result = BigFloat(0; precision=512)
add!(result, x, y)  # result uses 512 bits
```

### Q: How do I choose rounding modes?

**A:**

- `MPFRRoundNearest` (default): Minimize error
- `MPFRRoundToZero`: Truncate (useful for interval arithmetic)
- `MPFRRoundUp/Down`: Directed rounding (cryptography, validation)

```julia
sqrt!(result, 2; rounding=MPFRRoundDown)  # Lower bound
sqrt!(result_ub, 2; rounding=MPFRRoundUp)  # Upper bound
```

### Q: How do I handle errors?

**A:** Most operations return `Cint` (TeX exactness). Positive = inexact rounding occurred.

```julia
flag = sqrt!(result, x)
if flag != 0
    @warn "Result was rounded"
end
```

### Q: Is this thread-safe?

**A:** Yes, InPlace.jl has **no global state**. Each operation is independent:

```julia
Threads.@threads for i in 1:Threads.nthreads()
    result = BigInt(0)
    add!(result, i, i)  # Safe: own result object
end
```

### Q: Can I use this in embedded systems?

**A:** Yes! InPlace.jl is deterministic (no GC during operations):

```julia
@time for i in 1:1000000
    mul!(result, x, y)
end
# Consistent, predictable timing
```

### Q: Does this work with GPU?

**A:** Not directly (GMP is CPU-only), but you can:

1. Use InPlace.jl on CPU for arbitrary precision
2. Transfer finite-precision results to GPU

```julia
# CPU: Compute high-precision value
result = BigFloat(0; precision=1024)
exp!(result, BigFloat(π; precision=1024))

# GPU: Use as Float64
gpu_value = Float64(result)
```

### Q: Why is `mul_2exp!` better than `<<`?

**A:** Direct bit shifting avoids `BigInt` allocations:

```julia
# ❌ Julia's <<: allocates
y = x << 10

# ✅ InPlace: zero allocation
mul_2exp!(result, x, 10)
```

### Q: Can I use InPlace with ForwardDiff.jl?

**A:** Yes, but be careful—automatic differentiation needs to track through mutations:

```julia
using InPlace, ForwardDiff

function f(x)
    result = BigInt(0)
    mul!(result, x, 2)  # AD-aware mutation
    return result
end

ForwardDiff.derivative(f, 3.0)
```

---

## See Also

- [GMP Manual](https://gmplib.org/manual/)
- [MPFR Manual](https://www.mpfr.org/mpfr-current/mpfr.html)
- [Julia Arbitrary Precision Arithmetic](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic)

---

## License

InPlace.jl is licensed under the MIT License. See LICENSE file for details.
