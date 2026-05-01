using InPlace
using Test

function mpz_string(x::BigInt; base::Integer=10)
    n = InPlace.mpz_sizeinbase(x, base) + 2
    buf = Vector{UInt8}(undef, n)
    GC.@preserve buf begin
        ptr = InPlace.mpz_get_str(pointer(buf), base, x)
        @test ptr == pointer(buf)
        return unsafe_string(ptr)
    end
end

@testset "mpz assignment and conversion" begin
    x = BigInt(0)
    @test InPlace.mpz_set_str(x, "123456789012345678901234567890", 10) == 0
    @test x == big"123456789012345678901234567890"
    @test mpz_string(x) == "123456789012345678901234567890"

    @test InPlace.mpz_set_str(x, "-ff", 16) == 0
    @test x == -255
    @test mpz_string(x; base=16) == "-ff"
end

@testset "C integer range checks" begin
    x = BigInt(0)
    @test_throws InexactError InPlace.mpz_set_str(x, "10", Int64(typemax(Cint)) + 1)
    @test_throws InexactError InPlace.mpz_realloc2(x, -1)

    if Culong != UInt64
        @test_throws InexactError InPlace.mpz_set_ui(x, UInt64(typemax(Culong)) + 1)
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
    @test InPlace.mpz_mod_ui(rop, y, 97) == mod(y, UInt(97))
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
        InPlace.mpz_import(x, length(bytes), 1, 1, 1, 0, pointer(bytes))
    end
    @test x == 0x12345678

    out = Vector{UInt8}(undef, length(bytes))
    count = Ref{Csize_t}(0)
    GC.@preserve out begin
        ptr = InPlace.mpz_export(pointer(out), count, 1, 1, 1, 0, x)
        @test ptr == pointer(out)
    end
    @test count[] == length(bytes)
    @test out == bytes
end

@testset "mpz limbs" begin
    x = BigInt(0)
    InPlace.mpz_realloc2(x, 128)
    limbs = InPlace.mpz_limbs_write(x, 2)
    @test limbs != C_NULL
    unsafe_store!(limbs, InPlace.mp_limb_t(0x0000000000000001), 1)
    unsafe_store!(limbs, InPlace.mp_limb_t(0x0000000000000002), 2)
    InPlace.mpz_limbs_finish(x, 2)

    @test InPlace.mpz_size(x) == 2
    @test InPlace.mpz_getlimbn(x, 0) == InPlace.mp_limb_t(1)
    @test InPlace.mpz_getlimbn(x, 1) == InPlace.mp_limb_t(2)
    @test InPlace.mpz_limbs_read(x) != C_NULL
end

@testset "user-facing BigInt arithmetic" begin
    z = BigInt(0)
    same = z

    @test InPlace.set!(z, "12345678901234567890") === same
    @test z == big"12345678901234567890"

    @test InPlace.add!(z, z, 10) === same
    @test z == big"12345678901234567900"

    InPlace.sub!(z, z, 900)
    @test z == big"12345678901234567000"

    InPlace.mul!(z, z, -2)
    @test z == -big"24691357802469134000"

    InPlace.abs!(z, z)
    @test z == big"24691357802469134000"

    InPlace.div!(z, z, 10)
    @test z == big"2469135780246913400"

    InPlace.mod!(z, z, 97)
    @test z == mod(big"2469135780246913400", 97)

    InPlace.pow!(z, 12, 5)
    @test z == 12^5

    InPlace.gcd!(z, 48, 180)
    @test z == 12
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

@testset "user-facing BigFloat arithmetic" begin
    z = BigFloat(0; precision=256)
    same = z

    @test InPlace.set!(z, "1.25") === same
    @test z == BigFloat("1.25"; precision=256)

    @test InPlace.add!(z, z, 2) === same
    @test z == BigFloat("3.25"; precision=256)

    InPlace.mul!(z, z, 4)
    @test z == BigFloat(13; precision=256)

    InPlace.div!(z, z, 2)
    @test z == BigFloat("6.5"; precision=256)

    InPlace.sqrt!(z, 2)
    @test abs(z - sqrt(BigFloat(2; precision=256))) < big"1e-70"

    InPlace.fma!(z, 2, 3, 4)
    @test z == BigFloat(10; precision=256)

    InPlace.pi!(z)
    @test abs(z - BigFloat(pi; precision=256)) < big"1e-70"
end
