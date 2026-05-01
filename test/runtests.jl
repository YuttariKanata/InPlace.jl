using InPlace
using Test

function mpz_string(x::BigInt; base::Integer=10)
    0 <= base <= 62 ? base=Cint(base) : throw(ArgumentError("invalid integer string for base $base: $s"))
    n = InPlace.mpz_sizeinbase(x, base) + 2
    buf = Vector{UInt8}(undef, n)
    GC.@preserve buf begin
        ptr = InPlace.mpz_get_str(pointer(buf), base, x)
        @test ptr == pointer(buf)
        return unsafe_string(ptr)
    end
end

@testset "Low-level GMP/MPFR wrappers" begin

@testset "mpz assignment and conversion" begin
    x = BigInt(0)
    @test InPlace.set!(x, "123456789012345678901234567890";base=10) != 0
    @test x == big"123456789012345678901234567890"
    @test mpz_string(x) == "123456789012345678901234567890"

    @test InPlace.set!(x, "-ff", base=16) != 0
    @test x == -255
    @test mpz_string(x; base=16) == "-ff"
end

@testset "C integer range checks" begin
    x = BigInt(0)
    @test_throws ArgumentError InPlace.set!(x, "10"; base=63)
    @test_throws InexactError  InPlace.mpz_realloc2(x, -1)

    if Culong != UInt64
        @test_throws MethodError InPlace.mpz_set_ui(x, UInt128(10))
    end
end

@testset "mpz arithmetic and macro aliases" begin
    rop = BigInt(0)
    x = big"12345678901234567890"
    y = big"98765432109876543210"

    InPlace.mpz_add(rop, x, y)
    @test rop == x + y

    InPlace.mpz_mul(rop, x, y)
    @test rop == x * y

    InPlace.mpz_mod(rop, y, x)
    @test rop == mod(y, x)
    @test InPlace.mpz_mod_ui(rop, y, Culong(97)) == mod(y, UInt(97))
    @test rop == mod(y, 97)

    @test InPlace.mpz_sgn(BigInt(-3)) == -1
    @test InPlace.mpz_sgn(BigInt(0)) == 0
    @test InPlace.mpz_sgn(BigInt(3)) == 1
    @test InPlace.mpz_odd_p(BigInt(5)) != 0
    @test InPlace.mpz_even_p(BigInt(6)) != 0
    @test InPlace.mpz_kronecker(BigInt(5), BigInt(7)) == InPlace.mpz_jacobi(BigInt(5), BigInt(7))
end

@testset "mpz import/export" begin
    bytes = UInt8[0x12, 0x34, 0x56, 0x78]
    x = BigInt(0)
    GC.@preserve bytes begin
        InPlace.mpz_import(x, Csize_t(length(bytes)), Cint(1), Csize_t(1), Cint(1), Csize_t(0), pointer(bytes))
    end
    @test x == 0x12345678

    out = Vector{UInt8}(undef, length(bytes))
    count = Ref{Csize_t}(0)
    GC.@preserve out begin
        ptr = InPlace.mpz_export(pointer(out), count, Cint(1), Csize_t(1), Cint(1), Csize_t(0), x)
        @test ptr == pointer(out)
    end
    @test count[] == length(bytes)
    @test out == bytes
end

@testset "mpz limbs" begin
    x = BigInt(0)
    InPlace.mpz_realloc2(x, 128)
    limbs = InPlace.mpz_limbs_write(x, InPlace.mp_size_t(2))
    @test limbs != C_NULL
    unsafe_store!(limbs, InPlace.mp_limb_t(0x0000000000000001), 1)
    unsafe_store!(limbs, InPlace.mp_limb_t(0x0000000000000002), 2)
    InPlace.mpz_limbs_finish(x, InPlace.mp_size_t(2))

    @test InPlace.mpz_size(x) == 2
    @test InPlace.mpz_getlimbn(x, InPlace.mp_size_t(0)) == InPlace.mp_limb_t(1)
    @test InPlace.mpz_getlimbn(x, InPlace.mp_size_t(1)) == InPlace.mp_limb_t(2)
    @test InPlace.mpz_limbs_read(x) != C_NULL
end

@testset "mpfr assignment and conversion" begin
    x = BigFloat(0; precision=256)
    @test InPlace.mpfr_set_str(x, "1.25", 10, InPlace.MPFRRoundNearest) == 0
    @test x == BigFloat("1.25"; precision=256)
    @test InPlace.mpfr_get_d(x, InPlace.MPFRRoundNearest) == 1.25
    @test InPlace.mpfr_get_si(x, InPlace.MPFRRoundToZero) == 1

    z = BigInt(0)
    @test InPlace.mpfr_get_z(z, x, InPlace.MPFRRoundToZero) < 0
    @test z == 1

    exp = Ref{InPlace.mpfr_exp_t}(0)
    ptr = InPlace.mpfr_get_str(Ptr{UInt8}(C_NULL), exp, 10, 0, x, InPlace.MPFRRoundNearest)
    @test ptr != C_NULL
    s = unsafe_string(ptr)
    InPlace.mpfr_free_str(ptr)
    @test startswith(s, "125")
    @test exp[] == 1
end

@testset "mpfr arithmetic and elementary functions" begin
    x = BigFloat(2; precision=256)
    y = BigFloat(3; precision=256)
    z = BigFloat(0; precision=256)

    @test InPlace.mpfr_add(z, x, y, InPlace.MPFRRoundNearest) == 0
    @test z == BigFloat(5; precision=256)

    @test InPlace.mpfr_mul_ui(z, z, 4, InPlace.MPFRRoundNearest) == 0
    @test z == BigFloat(20; precision=256)

    InPlace.mpfr_sqrt(z, BigFloat(2; precision=256), InPlace.MPFRRoundNearest)
    @test abs(z - sqrt(BigFloat(2; precision=256))) < big"1e-70"

    InPlace.mpfr_const_pi(z, InPlace.MPFRRoundNearest)
    @test abs(z - BigFloat(pi; precision=256)) < big"1e-70"

    @test InPlace.mpfr_cmp_ui(BigFloat(10; precision=256), 10) == 0
    @test InPlace.mpfr_number_p(z) != 0
    @test InPlace.mpfr_nan_p(z) == 0

    InPlace.mpfr_min(z, BigFloat(2; precision=256), BigFloat(3; precision=256), InPlace.MPFRRoundNearest)
    @test z == 2
    InPlace.mpfr_max(z, BigFloat(2; precision=256), BigFloat(3; precision=256), InPlace.MPFRRoundNearest)
    @test z == 3
end

@testset "mpfr flags and range checks" begin
    InPlace.mpfr_clear_flags()
    @test InPlace.mpfr_inexflag_p() == 0
    InPlace.mpfr_set_divby0()
    @test InPlace.mpfr_divby0_p() != 0
    InPlace.mpfr_clear_divby0()
    @test InPlace.mpfr_divby0_p() == 0

    x = BigFloat(0; precision=256)
    @test_throws InexactError InPlace.mpfr_set_str(x, "1.0", Int64(typemax(Cint)) + 1, InPlace.MPFRRoundNearest)

    if Culong != UInt64
        @test_throws InexactError InPlace.mpfr_set_ui(x, UInt64(typemax(Culong)) + 1, InPlace.MPFRRoundNearest)
    end

    @test :mpfr_set_uj in InPlace.MPFR_MISSING
end

end # Low-level GMP/MPFR wrappers

##################################################################################################################################

@testset "BigInt API" begin

# ============================================================================
# BigInt: Assignment and Conversion Tests
# ============================================================================

@testset "BigInt: set! - String conversion" begin
    x = BigInt(0)
    
    # Decimal
    @test InPlace.set!(x, "123456789012345678901234567890") === x
    @test x == big"123456789012345678901234567890"
    
    # Hexadecimal
    InPlace.set!(x, "deadbeef", base=16)
    @test x == 0xdeadbeef
    
    # Binary
    InPlace.set!(x, "1010", base=2)
    @test x == 10
    
    # Negative
    InPlace.set!(x, "-999")
    @test x == -999
end

@testset "BigInt: set! - Integer conversion" begin
    x = BigInt(0)
    
    # From Int32
    InPlace.set!(x, Int32(42))
    @test x == 42
    
    # From Int64
    InPlace.set!(x, Int64(-999999999))
    @test x == -999999999
    
    # From UInt32
    InPlace.set!(x, UInt32(12345))
    @test x == 12345
    
    # From UInt64
    InPlace.set!(x, UInt64(18446744073709551615))
    @test x == 18446744073709551615
    
    # From BigInt
    y = BigInt(777)
    InPlace.set!(x, y)
    @test x == 777
end

@testset "BigInt: set! - Float conversion" begin
    x = BigInt(0)
    
    InPlace.set!(x, 3.7)
    @test x == 3
    
    InPlace.set!(x, -2.1)
    @test x == -2
end

@testset "BigInt: swap!" begin
    x = BigInt(100)
    y = BigInt(200)
    
    InPlace.swap!(x, y)
    @test x == 200
    @test y == 100
end

# ============================================================================
# BigInt: Addition Tests
# ============================================================================

@testset "BigInt: add! - BigInt + BigInt" begin
    result = BigInt(0)
    x = BigInt(100)
    y = BigInt(50)
    
    InPlace.add!(result, x, y)
    @test result == 150
end

@testset "BigInt: add! - BigInt + Native UnsignedInt" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    InPlace.add!(result, x, UInt32(42))
    @test result == 1042
    
    InPlace.add!(result, x, UInt64(999))
    @test result == 1999
end

@testset "BigInt: add! - BigInt + Native SignedInt" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    # Positive signed
    InPlace.add!(result, x, Int64(42))
    @test result == 1042
    
    # Negative signed
    InPlace.add!(result, x, Int64(-100))
    @test result == 900
end

@testset "BigInt: add! - Native + BigInt (commutative)" begin
    result = BigInt(0)
    x = BigInt(500)
    
    InPlace.add!(result, UInt64(250), x)
    @test result == 750
end

@testset "BigInt: add! - BigInt + Integer" begin
    result = BigInt(0)
    x = BigInt(100)
    y = Int16(50)
    
    InPlace.add!(result, x, y)
    @test result == 150
end

@testset "BigInt: add! - Integer + Integer" begin
    result = BigInt(0)
    
    InPlace.add!(result, 1000, 500)
    @test result == 1500
end

@testset "BigInt: operations use C scalar range for wide integer values" begin
    result = BigInt(0)
    x = BigInt(1000)

    small_signed = Int64(typemax(Clong))
    small_unsigned = UInt64(typemax(Culong))
    wide_unsigned = UInt64(typemax(Culong)) + UInt64(1)

    InPlace.set!(result, small_signed)
    @test result == small_signed

    InPlace.add!(result, x, small_signed)
    @test result == x + small_signed

    InPlace.sub!(result, x, small_signed)
    @test result == x - small_signed

    InPlace.mul!(result, x, small_signed)
    @test result == x * small_signed

    InPlace.add!(result, x, small_unsigned)
    @test result == x + small_unsigned

    InPlace.add!(result, x, wide_unsigned)
    @test result == x + wide_unsigned
end

# ============================================================================
# BigInt: Subtraction Tests
# ============================================================================

@testset "BigInt: sub! - BigInt - BigInt" begin
    result = BigInt(0)
    x = BigInt(1000)
    y = BigInt(300)
    
    InPlace.sub!(result, x, y)
    @test result == 700
end

@testset "BigInt: sub! - BigInt - Native UnsignedInt" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    InPlace.sub!(result, x, UInt64(100))
    @test result == 900
end

@testset "BigInt: sub! - BigInt - Native SignedInt" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    # Positive: normal subtraction
    InPlace.sub!(result, x, Int64(200))
    @test result == 800
    
    # Negative: becomes addition
    InPlace.sub!(result, x, Int64(-100))
    @test result == 1100
end

@testset "BigInt: sub! - Native UnsignedInt - BigInt" begin
    result = BigInt(0)
    y = BigInt(30)
    
    InPlace.sub!(result, UInt64(100), y)
    @test result == 70
end

@testset "BigInt: sub! - Native SignedInt - BigInt" begin
    result = BigInt(0)
    y = BigInt(50)
    
    # Positive: normal subtraction
    InPlace.sub!(result, Int64(200), y)
    @test result == 150
    
    # Negative: becomes negation of sum
    InPlace.sub!(result, Int64(-100), y)
    @test result == -150
end

@testset "BigInt: sub! - Large number" begin
    result = BigInt(0)
    x = big"999999999999999999999999999999"
    y = big"1"
    
    InPlace.sub!(result, x, y)
    @test result == big"999999999999999999999999999998"
end

# ============================================================================
# BigInt: Multiplication Tests
# ============================================================================

@testset "BigInt: mul! - BigInt * BigInt" begin
    result = BigInt(0)
    x = BigInt(12345)
    y = BigInt(67890)
    
    InPlace.mul!(result, x, y)
    @test result == 12345 * 67890
end

@testset "BigInt: mul! - BigInt * Native UnsignedInt" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    InPlace.mul!(result, x, UInt64(7))
    @test result == 7000
end

@testset "BigInt: mul! - BigInt * Native SignedInt" begin
    result = BigInt(0)
    x = BigInt(100)
    
    # Positive
    InPlace.mul!(result, x, Int64(5))
    @test result == 500
    
    # Negative
    InPlace.mul!(result, x, Int64(-3))
    @test result == -300
end

@testset "BigInt: mul! - Zero handling" begin
    result = BigInt(999)
    x = BigInt(555)
    
    InPlace.mul!(result, x, 0)
    @test result == 0
end

@testset "BigInt: mul! - Large multiplication" begin
    result = BigInt(0)
    x = big"123456789123456789"
    y = big"987654321987654321"
    
    InPlace.mul!(result, x, y)
    @test result == x * y
end

# ============================================================================
# BigInt: Division Tests
# ============================================================================

@testset "BigInt: div! - Truncated division" begin
    result = BigInt(0)
    x = BigInt(1000)
    y = BigInt(7)
    
    InPlace.div!(result, x, y; rounding=:trunc)
    @test result == 142
end

@testset "BigInt: div! - Floor division" begin
    result = BigInt(0)
    
    # Positive / Positive
    InPlace.div!(result, BigInt(1000), BigInt(7); rounding=:floor)
    @test result == 142
    
    # Negative / Positive
    InPlace.div!(result, BigInt(-1000), BigInt(7); rounding=:floor)
    @test result == -143  # Floors toward -∞
end

@testset "BigInt: div! - Ceiling division" begin
    result = BigInt(0)
    
    InPlace.div!(result, BigInt(1000), BigInt(7); rounding=:ceil)
    @test result == 143
end

@testset "BigInt: div! - Native UnsignedInt divisor" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    InPlace.div!(result, x, UInt64(7))
    @test result == 142
end

@testset "BigInt: div! - With remainder object" begin
    q = BigInt(0)
    r = BigInt(0)
    x = BigInt(1000)
    y = BigInt(7)
    
    InPlace.divrem!(q, r, x, y; rounding=:trunc)
    @test q == 142
    @test r == 6
    @test x == q * y + r
end

# ============================================================================
# BigInt: Remainder Tests
# ============================================================================

@testset "BigInt: rem! - Truncated remainder" begin
    result = BigInt(0)
    
    InPlace.rem!(result, BigInt(1000), BigInt(7); rounding=:trunc)
    @test result == 6
end

@testset "BigInt: rem! - Floor remainder" begin
    result = BigInt(0)
    
    InPlace.rem!(result, BigInt(-1000), BigInt(7); rounding=:floor)
    @test result == 1  # Different from truncated
end

@testset "BigInt: mod! - Euclidean modulo" begin
    result = BigInt(0)
    
    InPlace.mod!(result, BigInt(1000), BigInt(7))
    @test result == 6
    
    InPlace.mod!(result, BigInt(-1000), BigInt(7))
    @test result == 1  # Always non-negative
end

# ============================================================================
# BigInt: Compound Operations (addmul, submul, mul_2exp)
# ============================================================================

@testset "BigInt: addmul! - result += x * y" begin
    result = BigInt(100)
    x = BigInt(5)
    y = BigInt(10)
    
    InPlace.addmul!(result, x, y)
    @test result == 150  # 100 + 5*10
end

@testset "BigInt: addmul! - With native UnsignedInt" begin
    result = BigInt(100)
    x = BigInt(5)
    
    InPlace.addmul!(result, x, UInt64(20))
    @test result == 200  # 100 + 5*20
end

@testset "BigInt: addmul! - With native SignedInt (negative)" begin
    result = BigInt(100)
    x = BigInt(5)
    
    InPlace.addmul!(result, x, Int64(-10))
    @test result == 50  # 100 - 5*10
end

@testset "BigInt: submul! - result -= x * y" begin
    result = BigInt(100)
    x = BigInt(5)
    y = BigInt(10)
    
    InPlace.submul!(result, x, y)
    @test result == 50  # 100 - 5*10
end

@testset "BigInt: submul! - With native SignedInt (negative)" begin
    result = BigInt(100)
    x = BigInt(5)
    
    InPlace.submul!(result, x, Int64(-10))
    @test result == 150  # 100 + 5*10
end

@testset "BigInt: mul_2exp! - Bit shifting (multiply by 2^n)" begin
    result = BigInt(0)
    x = BigInt(7)
    
    InPlace.mul_2exp!(result, x, 4)
    @test result == 7 * 16  # 7 * 2^4
    @test result == 112
end

@testset "BigInt: div_2exp! - Bit shifting (divide by 2^n)" begin
    result = BigInt(0)
    x = BigInt(1000)
    
    InPlace.div_2exp!(result, x, 3; rounding=:trunc)
    @test result == 125  # 1000 / 2^3
end

# ============================================================================
# BigInt: Exponentiation Tests
# ============================================================================

@testset "BigInt: pow! - Small exponents" begin
    result = BigInt(0)
    
    InPlace.pow!(result, BigInt(2), UInt64(10))
    @test result == 1024
    
    InPlace.pow!(result, BigInt(3), UInt64(5))
    @test result == 243
end

@testset "BigInt: pow! - Large exponents" begin
    result = BigInt(0)
    
    InPlace.pow!(result, BigInt(2), UInt64(100))
    @test result == big"1267650600228229401496703205376"
end

@testset "BigInt: pow! - Native base and exponent" begin
    result = BigInt(0)
    
    InPlace.pow!(result, UInt64(2), UInt64(20))
    @test result == 1048576
end

@testset "BigInt: pow! - Zero exponent" begin
    result = BigInt(0)
    
    InPlace.pow!(result, 999, UInt64(0))
    @test result == 1
end

@testset "BigInt: powm! - Modular exponentiation" begin
    result = BigInt(0)
    base = BigInt(2)
    exp = BigInt(10)
    mod = BigInt(1000)
    
    InPlace.powm!(result, base, exp, mod)
    @test result == 24  # (2^10) mod 1000
end

@testset "BigInt: powm_sec! - Secure version" begin
    result = BigInt(0)
    base = BigInt(3)
    exp = BigInt(7)
    mod = BigInt(11)
    
    InPlace.powm_sec!(result, base, exp, mod)
    @test result == 9  # (3^7) mod 11
end

# ============================================================================
# BigInt: Root Tests
# ============================================================================

@testset "BigInt: sqrt! - Perfect squares" begin
    result = BigInt(0)
    
    InPlace.sqrt!(result, BigInt(16))
    @test result == 4
    
    InPlace.sqrt!(result, BigInt(10000))
    @test result == 100
end

@testset "BigInt: sqrt! - Non-perfect squares (floor)" begin
    result = BigInt(0)
    
    InPlace.sqrt!(result, BigInt(17))
    @test result == 4  # Floor(√17)
    
    InPlace.sqrt!(result, BigInt(100))
    @test result == 10
end

@testset "BigInt: sqrtrem! - Quotient and remainder" begin
    root = BigInt(0)
    rem = BigInt(0)
    x = BigInt(17)
    
    InPlace.sqrtrem!(root, rem, x)
    @test root == 4
    @test rem == 1
    @test 17 == root^2 + rem
end

@testset "BigInt: root! - N-th root" begin
    result = BigInt(0)
    
    InPlace.root!(result, BigInt(8), UInt64(3))
    @test result == 2  # ∛8
    
    InPlace.root!(result, BigInt(81), UInt64(4))
    @test result == 3  # ⁴√81
end

@testset "BigInt: rootrem! - N-th root with remainder" begin
    root = BigInt(0)
    rem = BigInt(0)
    
    InPlace.rootrem!(root, rem, BigInt(10), UInt64(3))
    @test root == 2
    @test rem == 2
    @test 10 == root^3 + rem
end

# ============================================================================
# BigInt: Unary Operations
# ============================================================================

@testset "BigInt: neg! - Negation" begin
    result = BigInt(0)
    
    InPlace.neg!(result, BigInt(42))
    @test result == -42
    
    InPlace.neg!(result, BigInt(-100))
    @test result == 100
end

@testset "BigInt: abs! - Absolute value" begin
    result = BigInt(0)
    
    InPlace.abs!(result, BigInt(42))
    @test result == 42
    
    InPlace.abs!(result, BigInt(-42))
    @test result == 42
end

# ============================================================================
# BigInt: Number Theoretic Functions
# ============================================================================

@testset "BigInt: gcd! - Greatest common divisor" begin
    result = BigInt(0)
    
    InPlace.gcd!(result, BigInt(48), BigInt(180))
    @test result == 12
    
    InPlace.gcd!(result, BigInt(100), BigInt(35))
    @test result == 5
end

@testset "BigInt: gcd! - With native UnsignedInt" begin
    result = BigInt(0)
    
    InPlace.gcd!(result, BigInt(48), UInt64(18))
    @test result == 6
end

@testset "BigInt: lcm! - Least common multiple" begin
    result = BigInt(0)
    
    InPlace.lcm!(result, BigInt(48), BigInt(180))
    @test result == 720
end

@testset "BigInt: gcdext! - Extended GCD" begin
    g = BigInt(0)
    s = BigInt(0)
    t = BigInt(0)
    
    InPlace.gcdext!(g, s, t, BigInt(35), BigInt(15))
    @test g == 5
    @test 35 * s + 15 * t == g
end

@testset "BigInt: invert! - Modular inverse" begin
    result = BigInt(0)
    
    # 3 * 5 ≡ 1 (mod 7)
    InPlace.invert!(result, BigInt(3), BigInt(7))
    @test result == 5
    @test (3 * result) % 7 == 1
end

@testset "BigInt: remove! - Remove factor" begin
    result = BigInt(0)
    
    # 72 = 2^3 * 3^2, remove factor 2
    count = InPlace.remove!(result, BigInt(72), BigInt(2))
    @test result == 9  # 72 / 2^3
    @test count == 3   # Removed 3 times
end

@testset "BigInt: fac_ui! - Factorial" begin
    result = BigInt(0)
    
    InPlace.fac_ui!(result, UInt64(5))
    @test result == 120
    
    InPlace.fac_ui!(result, UInt64(10))
    @test result == 3628800
end

@testset "BigInt: bin_ui! - Binomial coefficient" begin
    result = BigInt(0)
    
    InPlace.bin_ui!(result, BigInt(10), UInt64(3))
    @test result == 120  # C(10,3)
    
    InPlace.bin_ui!(result, BigInt(5), UInt64(2))
    @test result == 10   # C(5,2)
end

@testset "BigInt: bin_uiui! - Binomial (both args unsigned)" begin
    result = BigInt(0)
    
    InPlace.bin_uiui!(result, UInt64(10), UInt64(3))
    @test result == 120
end

@testset "BigInt: fib_ui! - Fibonacci number" begin
    result = BigInt(0)
    
    InPlace.fib_ui!(result, UInt64(10))
    @test result == 55
    
    InPlace.fib_ui!(result, UInt64(15))
    @test result == 610
end

@testset "BigInt: fib2_ui! - Fibonacci pair" begin
    f_n = BigInt(0)
    f_n1 = BigInt(0)
    
    InPlace.fib2_ui!(f_n, f_n1, UInt64(10))
    @test f_n == 55
    @test f_n1 == 34
end

@testset "BigInt: lucnum_ui! - Lucas number" begin
    result = BigInt(0)
    
    InPlace.lucnum_ui!(result, UInt64(5))
    @test result == 11
end

@testset "BigInt: lucnum2_ui! - Lucas number pair" begin
    l_n = BigInt(0)
    l_n1 = BigInt(0)
    
    InPlace.lucnum2_ui!(l_n, l_n1, UInt64(5))
    @test l_n == 11
    @test l_n1 == 7
end

# ============================================================================
# BigInt: Bitwise Operations
# ============================================================================

@testset "BigInt: and! - Bitwise AND" begin
    result = BigInt(0)
    
    InPlace.and!(result, BigInt(12), BigInt(10))
    @test result == 8  # 1100 AND 1010 = 1000
end

@testset "BigInt: or! - Bitwise OR" begin
    result = BigInt(0)
    
    InPlace.or!(result, BigInt(12), BigInt(10))
    @test result == 14  # 1100 OR 1010 = 1110
end

@testset "BigInt: xor! - Bitwise XOR" begin
    result = BigInt(0)
    
    InPlace.xor!(result, BigInt(12), BigInt(10))
    @test result == 6  # 1100 XOR 1010 = 0110
end

@testset "BigInt: com! - Bitwise complement" begin
    result = BigInt(0)
    
    InPlace.com!(result, BigInt(5))
    @test result == -6  # ~5 in two's complement
end

@testset "BigInt: setbit! - Set bit at position" begin
    result = BigInt(0)
    InPlace.set!(result, BigInt(5))  # 0101
    
    InPlace.setbit!(result, 1)  # Set bit 1
    @test result == 7  # 0111
end

@testset "BigInt: clrbit! - Clear bit at position" begin
    result = BigInt(0)
    InPlace.set!(result, BigInt(7))  # 0111
    
    InPlace.clrbit!(result, 1)  # Clear bit 1
    @test result == 5  # 0101
end

@testset "BigInt: combit! - Toggle bit at position" begin
    result = BigInt(0)
    InPlace.set!(result, BigInt(5))  # 0101
    
    InPlace.combit!(result, 1)  # Toggle bit 1
    @test result == 7  # 0111
    
    InPlace.combit!(result, 1)  # Toggle back
    @test result == 5  # 0101
end

end # BigInt API

@testset "BigFloat API" begin

# ============================================================================
# BigFloat: Assignment Tests
# ============================================================================

@testset "BigFloat: set! - String conversion" begin
    x = BigFloat(0; precision=256)
    
    InPlace.set!(x, "3.14159265358979323846")
    @test abs(x - BigFloat("3.14159265358979323846"; precision=256)) < big"1e-20"
end

@testset "BigFloat: set! - Integer conversion" begin
    x = BigFloat(0; precision=256)
    
    InPlace.set!(x, 42)
    @test x == 42
    
    InPlace.set!(x, Int64(-999))
    @test x == -999
end

@testset "BigFloat: set! - Float conversion" begin
    x = BigFloat(0; precision=256)
    
    InPlace.set!(x, Float64(3.14159))
    @test abs(x - 3.14159) < 1e-5
end

@testset "BigFloat: set! - BigInt conversion" begin
    x = BigFloat(0; precision=256)
    
    InPlace.set!(x, BigInt(1000000000000))
    @test x == 1000000000000
end

@testset "BigFloat: set! - Copy with rounding" begin
    x = BigFloat(π; precision=256)
    y = BigFloat(0; precision=128)
    
    InPlace.set!(y, x; rounding=InPlace.MPFRRoundNearest)
    @test abs(y - π) < 1e-30
end

@testset "BigFloat: swap!" begin
    x = BigFloat(2.5; precision=256)
    y = BigFloat(3.5; precision=256)
    
    InPlace.swap!(x, y)
    @test x == 3.5
    @test y == 2.5
end

# ============================================================================
# BigFloat: Arithmetic Tests
# ============================================================================

@testset "BigFloat: add! - BigFloat + BigFloat" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(2.5; precision=256)
    y = BigFloat(1.5; precision=256)
    
    InPlace.add!(result, x, y)
    @test result == 4.0
end

@testset "BigFloat: add! - BigFloat + native UnsignedInt" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(2.5; precision=256)
    
    InPlace.add!(result, x, UInt64(1))
    @test result == 3.5
end

@testset "BigFloat: add! - BigFloat + native SignedInt" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(2.5; precision=256)
    
    InPlace.add!(result, x, Int64(-1))
    @test result == 1.5
end

@testset "BigFloat: sub! - BigFloat - BigFloat" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(5.0; precision=256)
    y = BigFloat(2.0; precision=256)
    
    InPlace.sub!(result, x, y)
    @test result == 3.0
end

@testset "BigFloat: mul! - BigFloat * BigFloat" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(2.5; precision=256)
    y = BigFloat(4.0; precision=256)
    
    InPlace.mul!(result, x, y)
    @test result == 10.0
end

@testset "BigFloat: div! - BigFloat / BigFloat" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(10.0; precision=256)
    y = BigFloat(2.5; precision=256)
    
    InPlace.div!(result, x, y)
    @test result == 4.0
end

# ============================================================================
# BigFloat: Unary Operations
# ============================================================================

@testset "BigFloat: sqrt! - Square root" begin
    result = BigFloat(0; precision=256)
    
    InPlace.sqrt!(result, 2)
    @test abs(result - sqrt(BigFloat(2; precision=256))) < big"1e-70"
end

@testset "BigFloat: sqrt! - Perfect square" begin
    result = BigFloat(0; precision=256)
    
    InPlace.sqrt!(result, 16)
    @test result == 4.0
end

@testset "BigFloat: square! - x^2" begin
    result = BigFloat(0; precision=256)
    
    InPlace.square!(result, 3)
    @test result == 9.0
end

@testset "BigFloat: cbrt! - Cube root" begin
    result = BigFloat(0; precision=256)
    
    InPlace.cbrt!(result, 8)
    @test result == 2.0
end

@testset "BigFloat: rec_sqrt! - 1/√x" begin
    result = BigFloat(0; precision=256)
    
    InPlace.rec_sqrt!(result, 4)
    @test result == 0.5
end

@testset "BigFloat: neg! - Negation" begin
    result = BigFloat(0; precision=256)
    
    InPlace.neg!(result, 42)
    @test result == -42
end

@testset "BigFloat: abs! - Absolute value" begin
    result = BigFloat(0; precision=256)
    
    InPlace.abs!(result, -3.14)
    @test result == 3.14
end

# ============================================================================
# BigFloat: Exponential Functions
# ============================================================================

@testset "BigFloat: exp! - e^x" begin
    result = BigFloat(0; precision=256)
    
    InPlace.exp!(result, 1)
    @test abs(result - exp(BigFloat(1; precision=256))) < big"1e-70"
end

@testset "BigFloat: exp2! - 2^x" begin
    result = BigFloat(0; precision=256)
    
    InPlace.exp2!(result, 3)
    @test result == 8.0
end

@testset "BigFloat: exp10! - 10^x" begin
    result = BigFloat(0; precision=256)
    
    InPlace.exp10!(result, 2)
    @test result == 100.0
end

@testset "BigFloat: expm1! - e^x - 1" begin
    result = BigFloat(0; precision=256)
    
    InPlace.expm1!(result, 1)
    @test abs(result - (exp(BigFloat(1; precision=256)) - 1)) < big"1e-70"
end

# ============================================================================
# BigFloat: Logarithmic Functions
# ============================================================================

@testset "BigFloat: log! - Natural logarithm" begin
    result = BigFloat(0; precision=256)
    
    InPlace.log!(result, exp(BigFloat(1; precision=256)))
    @test abs(result - 1) < big"1e-70"
end

@testset "BigFloat: log2! - Base 2 logarithm" begin
    result = BigFloat(0; precision=256)
    
    InPlace.log2!(result, 8)
    @test result == 3.0
end

@testset "BigFloat: log10! - Base 10 logarithm" begin
    result = BigFloat(0; precision=256)
    
    InPlace.log10!(result, 100)
    @test result == 2.0
end

@testset "BigFloat: log1p! - log(1+x)" begin
    result = BigFloat(0; precision=256)
    
    InPlace.log1p!(result, BigFloat(0; precision=256))
    @test result == 0.0
end

# ============================================================================
# BigFloat: Trigonometric Functions
# ============================================================================

@testset "BigFloat: sin! - Sine" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(0; precision=256)
    InPlace.pi!(x)
    
    InPlace.sin!(result, x)
    @test abs(result) < big"1e-70"  # sin(π) ≈ 0
end

@testset "BigFloat: cos! - Cosine" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(0; precision=256)
    InPlace.pi!(x)
    
    InPlace.cos!(result, x)
    @test abs(result + 1) < big"1e-70"  # cos(π) = -1
end

@testset "BigFloat: tan! - Tangent" begin
    result = BigFloat(0; precision=256)
    
    InPlace.tan!(result, 0)
    @test result == 0.0
end

@testset "BigFloat: sin_cos! - Simultaneous sin and cos" begin
    s = BigFloat(0; precision=256)
    c = BigFloat(0; precision=256)
    x = BigFloat(0; precision=256)
    InPlace.pi!(x)
    InPlace.div!(x, x, 4)  # π/4
    
    InPlace.sin_cos!(s, c, x)
    @test abs(s - c) < big"1e-70"  # sin(π/4) = cos(π/4)
end

# ============================================================================
# BigFloat: Special Functions
# ============================================================================

@testset "BigFloat: gamma! - Gamma function" begin
    result = BigFloat(0; precision=256)
    
    InPlace.gamma!(result, 5)
    @test abs(result - 24) < big"1e-70"  # Γ(5) = 4! = 24
end

@testset "BigFloat: zeta! - Riemann zeta" begin
    result = BigFloat(0; precision=256)
    
    InPlace.zeta!(result, 2)
    @test abs(result - (BigFloat(π; precision=256)^2 / 6)) < big"1e-70"  # ζ(2) = π²/6
end

@testset "BigFloat: erf! - Error function" begin
    result = BigFloat(0; precision=256)
    
    InPlace.erf!(result, 0)
    @test result == 0.0
end

@testset "BigFloat: erfc! - Complementary error function" begin
    result = BigFloat(0; precision=256)
    
    InPlace.erfc!(result, 0)
    @test result == 1.0
end

# ============================================================================
# BigFloat: Advanced Operations
# ============================================================================

@testset "BigFloat: fma! - Fused multiply-add (x*y+z)" begin
    result = BigFloat(0; precision=256)
    
    InPlace.fma!(result, 2, 3, 4)
    @test result == 10  # 2*3 + 4
end

@testset "BigFloat: fms! - Fused multiply-subtract (x*y-z)" begin
    result = BigFloat(0; precision=256)
    
    InPlace.fms!(result, 2, 3, 4)
    @test result == 2  # 2*3 - 4
end

@testset "BigFloat: hypot! - Euclidean distance √(x²+y²)" begin
    result = BigFloat(0; precision=256)
    
    InPlace.hypot!(result, 3, 4)
    @test result == 5.0
end

@testset "BigFloat: dim! - Positive difference max(x-y,0)" begin
    result = BigFloat(0; precision=256)
    
    InPlace.dim!(result, 5, 2)
    @test result == 3.0
    
    InPlace.dim!(result, 2, 5)
    @test result == 0.0
end

@testset "BigFloat: min! - Minimum" begin
    result = BigFloat(0; precision=256)
    
    InPlace.min!(result, 2, 3)
    @test result == 2
    
    InPlace.min!(result, 5, 1)
    @test result == 1
end

@testset "BigFloat: max! - Maximum" begin
    result = BigFloat(0; precision=256)
    
    InPlace.max!(result, 2, 3)
    @test result == 3
    
    InPlace.max!(result, 5, 1)
    @test result == 5
end

# ============================================================================
# BigFloat: Rounding Functions
# ============================================================================

@testset "BigFloat: ceil! - Round up" begin
    result = BigFloat(0; precision=256)
    
    InPlace.ceil!(result, 2.3)
    @test result == 3.0
    
    InPlace.ceil!(result, -2.3)
    @test result == -2.0
end

@testset "BigFloat: floor! - Round down" begin
    result = BigFloat(0; precision=256)
    
    InPlace.floor!(result, 2.7)
    @test result == 2.0
    
    InPlace.floor!(result, -2.7)
    @test result == -3.0
end

@testset "BigFloat: round! - Banker's rounding" begin
    result = BigFloat(0; precision=256)
    
    InPlace.round!(result, 2.5)
    @test result == 2.0 || result == 3.0  # Banker's rounding (round to even)
end

@testset "BigFloat: trunc! - Truncate toward zero" begin
    result = BigFloat(0; precision=256)
    
    InPlace.trunc!(result, 2.7)
    @test result == 2.0
    
    InPlace.trunc!(result, -2.7)
    @test result == -2.0
end

@testset "BigFloat: frac! - Fractional part" begin
    result = BigFloat(0; precision=256)
    
    InPlace.frac!(result, BigFloat("2.7"; precision=256))
    @test abs(result - BigFloat("0.7"; precision=256)) < big"1e-70"
end

@testset "BigFloat: modf! - Split integer and fractional parts" begin
    int_part = BigFloat(0; precision=256)
    frac_part = BigFloat(0; precision=256)
    
    InPlace.modf!(int_part, frac_part, BigFloat("2.7"; precision=256))
    @test int_part == 2.0
    @test abs(frac_part - BigFloat("0.7"; precision=256)) < big"1e-70"
end

# ============================================================================
# BigFloat: Mathematical Constants
# ============================================================================

@testset "BigFloat: pi! - Compute π" begin
    result = BigFloat(0; precision=256)
    
    InPlace.pi!(result)
    @test abs(result - BigFloat(π; precision=256)) < big"1e-70"
end

@testset "BigFloat: euler! - Euler-Mascheroni constant γ" begin
    result = BigFloat(0; precision=256)
    
    InPlace.euler!(result)
    @test result > 0.5  # γ ≈ 0.5772...
    @test result < 0.6
end

@testset "BigFloat: catalan! - Catalan's constant" begin
    result = BigFloat(0; precision=256)
    
    InPlace.catalan!(result)
    @test result > 0.9  # G ≈ 0.9159...
    @test result < 1.0
end

end # BigFloat API

@testset "Integration and error behavior" begin

# ============================================================================
# Rounding Mode Tests
# ============================================================================

@testset "BigFloat: Rounding modes" begin
    x = BigFloat(√2; precision=64)
    r_nearest = BigFloat(0; precision=32)
    r_down = BigFloat(0; precision=32)
    r_up = BigFloat(0; precision=32)
    
    InPlace.set!(r_nearest, x; rounding=InPlace.MPFRRoundNearest)
    InPlace.set!(r_down, x; rounding=InPlace.MPFRRoundDown)
    InPlace.set!(r_up, x; rounding=InPlace.MPFRRoundUp)
    
    @test r_down <= x
    @test r_up >= x
    @test r_down <= r_nearest <= r_up
end

# ============================================================================
# Large Number Tests
# ============================================================================

@testset "BigInt: Very large numbers" begin
    result = BigInt(0)
    x = big"999999999999999999999999999999999999999999999"
    y = big"111111111111111111111111111111111111111111111"
    
    InPlace.add!(result, x, y)
    @test result == x + y
    
    InPlace.mul!(result, x, y)
    @test result == x * y
end

@testset "BigFloat: High precision" begin
    result = BigFloat(0; precision=2048)
    
    InPlace.pi!(result)
    # Check that precision is actually being used
    @test precision(result) == 2048
end

# ============================================================================
# Error Handling Tests
# ============================================================================

@testset "BigFloat: Division by zero" begin
    result = BigFloat(0; precision=256)
    x = BigFloat(1; precision=256)
    y = BigFloat(0; precision=256)
    
    # MPFR should set infinity
    InPlace.div!(result, x, y)
    @test isinf(result)
end

@testset "BigInt: Modular inverse failure" begin
    result = BigInt(0)
    
    # No modular inverse for non-coprime numbers
    @test_throws ErrorException InPlace.invert!(result, BigInt(4), BigInt(6))
end

end # Integration and error behavior

println("All tests completed successfully!")
