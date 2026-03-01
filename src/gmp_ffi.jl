# src/gmp_ffi.jl
# Low-level GMP bindings using GMP_jll

using GMP_jll
using Base.GMP.MPZ: mpz_t, bitcnt_t
using Base.GMP: CulongMax, ClongMax, CdoubleMax

# --- Type Aliases ---
if Clong == Int32
    const LargeInt = Union{Int64, Int128}
    const LargeUInt = Union{UInt64, UInt128}
else
    const LargeInt = Int128
    const LargeUInt = UInt128
end

# Helper to construct ccall targets
gmpz(op::Symbol) = (Symbol(:__g, op), libgmp)

# Helper to map C types to Julia types for function signatures
function cnv(op::Symbol)::Symbol
    if op === :mpz_t
        return :BigInt
    elseif op === :bitcnt_t
        return :bitcnt_t
    elseif op === :Culong
        return :CulongMax
    elseif op === :Clong
        return :ClongMax
    elseif op === :Cdouble
        return :CdoubleMax
    else
        return Symbol(op, :Max)
    end
end

# --- Arithmetic: addmul, submul ---
for (fname, gmpname, ytype) in [
    (:addmul!, :mpz_addmul,    :mpz_t ),
    (:addmul!, :mpz_addmul_ui, :Culong),
    (:submul!, :mpz_submul,    :mpz_t ),
    (:submul!, :mpz_submul_ui, :Culong),]

    @eval begin
        function $fname(z::BigInt, x::BigInt, y::$(cnv(ytype)))::Nothing
            ccall($(gmpz(gmpname)), Cvoid, (mpz_t, mpz_t, $ytype), z, x, y)
            return nothing
        end
    end
end

# --- Number Theory: Division & Congruence ---
for (fname, gmpname, ytype) in [
    (:isdivisible,      :mpz_divisible_p,      :mpz_t   ),
    (:isdivisible,      :mpz_divisible_ui_p,   :Culong  ),
    (:isdivisible_2exp, :mpz_divisible_2exp_p, :bitcnt_t),]

    @eval begin
        function $fname(n::BigInt, d::$(cnv(ytype)))::Bool
            return !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t, $ytype), n, d))
        end
    end
end

for (fname, gmpname, ytype) in [
    (:divexact!, :mpz_divexact,    :mpz_t ),
    (:divexact!, :mpz_divexact_ui, :Culong),]

    @eval begin
        function $fname(z::BigInt, n::BigInt, d::$(cnv(ytype)))::BigInt
            ccall($(gmpz(gmpname)), Cvoid, (mpz_t, mpz_t, $ytype), z, n, d)
            return z
        end
    end
end

"""
    mod_ui(a::BigInt, b::CulongMax)::Culong

Return `a % b`. Wrapper for `__gmpz_fdiv_ui`.
Note: This returns a `Culong`, not a `BigInt`.
"""
function mod_ui(a::BigInt, b::CulongMax)::Culong
    return ccall((:__gmpz_fdiv_ui, libgmp), Culong, (mpz_t, Culong), a, b)
end

for (fname, gmpname, ytype1, ytype2) in [
    (:iscongruent,      :mpz_congruent_p,      :mpz_t,  :mpz_t   ),
    (:iscongruent,      :mpz_congruent_ui_p,   :Culong, :Culong  ),
    (:iscongruent_2exp, :mpz_congruent_2exp_p, :mpz_t,  :bitcnt_t),]

    @eval begin
        function $fname(n::BigInt, c::$(cnv(ytype1)), d::$(cnv(ytype2)))::Bool
            return !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t, $ytype1, $ytype2), n, c, d))
        end
    end
end

# --- Powers & Roots ---
function powm!(rop::BigInt, base::BigInt, exp::CulongMax, mod::BigInt)::BigInt
    ccall((:__gmpz_powm_ui, libgmp), Cvoid, (mpz_t, mpz_t, Culong, mpz_t), rop, base, exp, mod)
    return rop
end

function pow_ui!(rop::BigInt, base::CulongMax, exp::CulongMax)::BigInt
    ccall((:__gmpz_ui_pow_ui, libgmp), Cvoid, (mpz_t, Culong, Culong), rop, base, exp)
    return rop
end

function iroot!(z::BigInt, x::BigInt, n::CulongMax)::Bool
    # Returns true if the root is exact
    return !iszero(ccall((:__gmpz_root, libgmp), Cint, (mpz_t, mpz_t, Culong), z, x, n))
end

function rootrem!(root::BigInt, rem::BigInt, x::BigInt, n::CulongMax)::Nothing
    ccall((:__gmpz_rootrem, libgmp), Cvoid, (mpz_t, mpz_t, mpz_t, Culong), root, rem, x, n)
    return nothing
end

function sqrtrem!(root::BigInt, rem::BigInt, x::BigInt)::Nothing
    ccall((:__gmpz_sqrtrem, libgmp), Cvoid, (mpz_t, mpz_t, mpz_t), root, rem, x)
    return nothing
end

for (fname, gmpname) in [
    (:isperfectpower,  :mpz_perfect_power_p ),
    (:isperfectsquare, :mpz_perfect_square_p),]

    @eval begin
        function $fname(n::BigInt)::Bool
            return !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t,), n))
        end
    end
end

# --- GCD & LCM ---
function gcd(op1::BigInt, op2::CulongMax)::Culong
    # Note: Returning Culong directly as per GMP spec for _ui variants
    return ccall((:__gmpz_gcd_ui, libgmp), Culong, (Ptr{Cvoid}, mpz_t, Culong), C_NULL, op1, op2)
end

function lcm!(rop::BigInt, op1::BigInt, op2::CulongMax)::BigInt
    ccall((:__gmpz_lcm_ui, libgmp), Cvoid, (mpz_t, mpz_t, Culong), rop, op1, op2)
    return rop
end

# --- Number Theory Special Functions ---
for (fname, gmpname, ytype1, ytype2) in [
    (:jacobi,      :mpz_jacobi,       :mpz_t,  :mpz_t ),
    (:legendre,    :mpz_legendre,     :mpz_t,  :mpz_t ),
    (:kronecker,   :mpz_kronecker,    :mpz_t,  :mpz_t ),
    (:kronecker,   :mpz_kronecker_si, :mpz_t,  :Clong ),
    (:kronecker,   :mpz_kronecker_ui, :mpz_t,  :Culong),
    (:kronecker,   :mpz_si_kronecker, :Clong,  :mpz_t ),
    (:kronecker,   :mpz_ui_kronecker, :Culong, :mpz_t ),]

    @eval begin
        function $fname(n::$(cnv(ytype1)), k::$(cnv(ytype2)))::Cint
            return ccall($(gmpz(gmpname)), Cint, ($ytype1, $ytype2), n, k)
        end
    end
end

function removefactor!(rop::BigInt, op::BigInt, f::BigInt)::bitcnt_t
    return ccall((:__gmpz_remove, libgmp), bitcnt_t, (mpz_t, mpz_t, mpz_t), rop, op, f)
end

# --- Combinatorics ---
function fac2!(rop::BigInt, n::CulongMax)::BigInt
    ccall((:__gmpz_2fac_ui, libgmp), Cvoid, (mpz_t, Culong), rop, n)
    return rop
end

function facm!(rop::BigInt, n::CulongMax, m::CulongMax)::BigInt
    ccall((:__gmpz_mfac_uiui, libgmp), Cvoid, (mpz_t, Culong, Culong), rop, n, m)
    return rop
end

function primorial!(rop::BigInt, n::CulongMax)::BigInt
    ccall((:__gmpz_primorial_ui, libgmp), Cvoid, (mpz_t, Culong), rop, n)
    return rop
end

function binomial!(rop::BigInt, n::CulongMax, k::CulongMax)::BigInt
    ccall((:__gmpz_bin_uiui, libgmp), Cvoid, (mpz_t, Culong, Culong), rop, n, k)
    return rop
end

# --- Fibonacci & Lucas ---
for (fname, gmpname) in [(:fibonacci!, :mpz_fib_ui), (:lucas!, :mpz_lucnum_ui)]
    @eval begin
        function $fname(rop::BigInt, n::CulongMax)::BigInt
            ccall($(gmpz(gmpname)), Cvoid, (mpz_t, Culong), rop, n)
            return rop
        end
    end
end

for (fname, gmpname) in [
    (:fibonacci2!, :mpz_fib2_ui   ),
    (:lucas2!,     :mpz_lucnum2_ui),]
    @eval begin
        """
            $($fname)(rop::BigInt, ropsub1::BigInt, n::CulongMax)

        Calculate both `F_n` and `F_{n-1}` (or `L_n` and `L_{n-1}`) and store them in `rop` and `ropsub1`.
        """
        function $fname(rop::BigInt, ropsub1::BigInt, n::CulongMax)::Tuple{BigInt, BigInt}
            ccall($(gmpz(gmpname)), Cvoid, (mpz_t, mpz_t, Culong), rop, ropsub1, n)
            return (rop, ropsub1)
        end
    end
end

# --- Comparison & Misc ---
for (gmpname, ytype) in [
    (:mpz_cmpabs,    :mpz_t  ),
    (:mpz_cmpabs_d,  :Cdouble),
    (:mpz_cmpabs_ui, :Culong ),]
    
    @eval begin
        function cmpabs(op1::BigInt, op2::$(cnv(ytype)))::Cint
            return ccall($(gmpz(gmpname)), Cint, (mpz_t, $ytype), op1, op2)
        end
    end
end

function hamdist(op1::BigInt, op2::BigInt)::bitcnt_t
    return ccall((:__gmpz_hamdist, libgmp), bitcnt_t, (mpz_t, mpz_t), op1, op2)
end

function clrbit!(x::BigInt, a::bitcnt_t)::BigInt
    ccall((:__gmpz_clrbit, libgmp), Cvoid, (mpz_t, bitcnt_t), x, a)
    return x
end

function combit!(x::BigInt, a::bitcnt_t)::BigInt
    ccall((:__gmpz_combit, libgmp), Cvoid, (mpz_t, bitcnt_t), x, a)
    return x
end

for (gmpname, ytype) in [
    (:mpz_fits_ulong_p,  :Culong ),
    (:mpz_fits_slong_p,  :Clong  ),
    (:mpz_fits_uint_p,   :Cuint  ),
    (:mpz_fits_sint_p,   :Cint   ),
    (:mpz_fits_ushort_p, :Cushort),
    (:mpz_fits_sshort_p, :Cshort ),]

    @eval begin
        function isfit(op::BigInt, ::Type{$ytype})::Bool
            return !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t,), op))
        end
    end 
end