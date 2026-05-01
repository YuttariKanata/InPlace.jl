# User-facing mutating arithmetic helpers with Two-Tier Dispatch strategy
# Tier 1: Compile-time type dispatch (Static) via `@inline where {T <: <dynamic-native-type>}`
# Tier 2: Run-time value checks for fallback methods
# 
# API Design: Unified function names with type-driven dispatch
# (e.g., mul_2exp! handles both unsigned and signed exponents internally)

using Base.GMP: CulongMax, ClongMax

export set!,
    swap!,
    add!,
    sub!,
    mul!,
    div!,
    rem!,
    mod!,
    pow!,
    sqrt!,
    root!,
    neg!,
    abs!,
    gcd!,
    lcm!,
    and!,
    or!,
    xor!,
    com!,
    square!,
    fma!,
    fms!,
    sin!,
    cos!,
    tan!,
    acos!,
    asin!,
    atan!,
    atan2!,
    sinh!,
    cosh!,
    tanh!,
    asinh!,
    acosh!,
    atanh!,
    exp!,
    exp2!,
    exp10!,
    expm1!,
    log!,
    log2!,
    log10!,
    log1p!,
    min!,
    max!,
    pi!,
    euler!,
    catalan!,
    # Compound operations
    addmul!,
    submul!,
    mul_2exp!,
    div_2exp!,
    # Root extraction
    rootrem!,
    sqrtrem!,
    # Division with quotient/remainder
    divexact!,
    divrem!,
    divrem_ui!,
    # Exponentiation
    powm!,
    powm_sec!,
    # Number theoretic
    nextprime!,
    gcd_ui!,
    gcdext!,
    lcm_ui!,
    invert!,
    remove!,
    # Combinatorics
    fac_ui!,
    bin_ui!,
    bin_uiui!,
    fib_ui!,
    fib2_ui!,
    lucnum_ui!,
    lucnum2_ui!,
    # Bit manipulation
    setbit!,
    clrbit!,
    combit!,
    # MPFR specific (user-friendly unified names)
    dim!,
    rec_sqrt!,
    cbrt!,
    rootn!,
    sec!,
    csc!,
    cot!,
    sin_cos!,
    sinh_cosh!,
    sech!,
    csch!,
    coth!,
    gamma!,
    gamma_inc!,
    lngamma!,
    lgamma!,
    digamma!,
    beta!,
    zeta!,
    erf!,
    erfc!,
    j0!,
    j1!,
    jn!,
    y0!,
    y1!,
    yn!,
    agm!,
    ai!,
    eint!,
    li2!,
    rint!,
    ceil!,
    floor!,
    round!,
    trunc!,
    frac!,
    modf!,
    fmodquo!,
    remquo!,
    hypot!,
    reldiff!,
    nexttoward!,
    copysign!,
    nextabove!,
    nextbelow!,
    get_exp!,
    set_exp!,
    setsign!

# ============================================================================
# Type aliases and utilities (Tier 1/2 infrastructure)
# ============================================================================

const NativeSignedInt = ClongMax
const NativeUnsignedInt = CulongMax
const NativeInt = Union{NativeSignedInt, NativeUnsignedInt}

const FloatLike = Union{AbstractFloat, Integer, BigInt, BigFloat}

@inline _rounding(r::MPFRRoundingMode) = r
@inline _rounding(r::RoundingMode) = convert(MPFRRoundingMode, r)

@inline _as_bigint(x::BigInt) = x
@inline _as_bigint(x::Integer) = BigInt(x)

@inline _fits_clong(x::Integer) = typemin(Clong) <= x <= typemax(Clong)
@inline _fits_culong(x::Integer) = 0 <= x <= typemax(Culong)
@inline _as_clong(x::Integer) = Clong(x)
@inline _as_culong(x::Integer) = Culong(x)

@inline _as_bigfloat(x::BigFloat) = x
@inline _as_bigfloat(x::FloatLike; precision::Integer=precision(BigFloat)) =
    BigFloat(x; precision=precision)

@inline function _checked_nonnegative(::Val{name}, x::Integer) where {name}
    x < 0 && throw(ArgumentError("$(name) requires a non-negative integer argument"))
    return x
end

@inline function _with_bigfloat_arg(f, rop::BigFloat, x; rounding=MPFRRoundNearest)
    return f(rop, _as_bigfloat(x; precision=precision(rop)), _rounding(rounding))
end

# ============================================================================
# BigInt: Basic Arithmetic (Tier 1 + Tier 2 Hybrid)
# ============================================================================

@inline function set!(rop::BigInt, x::BigInt)
    mpz_set(rop, x)
    return rop
end

@inline function set!(rop::BigInt, x::T) where {T <: NativeSignedInt}
    mpz_set_si(rop, x)
    return rop
end

@inline function set!(rop::BigInt, x::T) where {T <: NativeUnsignedInt}
    mpz_set_ui(rop, x)
    return rop
end

@inline function set!(rop::BigInt, x::Integer)
    if x isa BigInt
        mpz_set(rop, x)
    elseif x isa NativeSignedInt
        mpz_set_si(rop, x)
    elseif x isa NativeUnsignedInt
        mpz_set_ui(rop, x)
    elseif _fits_clong(x)
        mpz_set_si(rop, _as_clong(x))
    elseif _fits_culong(x)
        mpz_set_ui(rop, _as_culong(x))
    else
        mpz_set(rop, _as_bigint(x))
    end
    return rop
end

function set!(rop::BigInt, s::AbstractString; base::Integer=0)
    0 <= base <= 62 ? base=Cint(base) : throw(ArgumentError("invalid integer string for base $base: $s"))
    rc = mpz_set_str(rop, s, base)
    rc == 0 || throw(ArgumentError("invalid integer string for base $base: $s"))
    return rop
end

@inline function set!(rop::BigInt, x::AbstractFloat)
    mpz_set_d(rop, x)
    return rop
end

@inline function swap!(x::BigInt, y::BigInt)
    mpz_swap(x, y)
    return x
end

# ============================================================================
# BigInt: add! (Unified multi-dispatch)
# ============================================================================

@inline function add!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_add(rop, x, y)
    return rop
end

@inline function add!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_add_ui(rop, x, y)
    return rop
end

@inline function add!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    if y < 0
        mpz_sub_ui(rop, x, Culong(-y))
    else
        mpz_add_ui(rop, x, Culong(y))
    end
    return rop
end

@inline function add!(rop::BigInt, x::T, y::BigInt) where {T <: NativeInt}
    return add!(rop, y, x)
end

@inline function add!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeInt
        add!(rop, x, y)
    elseif _fits_clong(y)
        add!(rop, x, _as_clong(y))
    elseif _fits_culong(y)
        add!(rop, x, _as_culong(y))
    else
        add!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function add!(rop::BigInt, x::Integer, y::BigInt)
    if x isa NativeInt
        add!(rop, y, x)
    elseif _fits_clong(x)
        add!(rop, y, _as_clong(x))
    elseif _fits_culong(x)
        add!(rop, y, _as_culong(x))
    else
        add!(rop, _as_bigint(x), y)
    end
    return rop
end

@inline function add!(rop::BigInt, x::T, y::U) where {T <: NativeInt, U <: NativeInt}
    return set!(rop, BigInt(x) + BigInt(y))
end

@inline function add!(rop::BigInt, x::Integer, y::Integer)
    return set!(rop, _as_bigint(x) + _as_bigint(y))
end

# ============================================================================
# BigInt: sub! (Unified multi-dispatch with mpz_ui_sub optimization)
# ============================================================================

@inline function sub!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_sub(rop, x, y)
    return rop
end

@inline function sub!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_sub_ui(rop, x, y)
    return rop
end

@inline function sub!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    if y < 0
        mpz_add_ui(rop, x, Culong(-y))
    else
        mpz_sub_ui(rop, x, Culong(y))
    end
    return rop
end

@inline function sub!(rop::BigInt, x::T, y::BigInt) where {T <: NativeUnsignedInt}
    mpz_ui_sub(rop, x, y)
    return rop
end

@inline function sub!(rop::BigInt, x::T, y::BigInt) where {T <: NativeSignedInt}
    if x < 0
        mpz_add_ui(rop, y, Culong(-x))
        mpz_neg(rop, rop)
    else
        mpz_ui_sub(rop, Culong(x), y)
    end
    return rop
end

@inline function sub!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeInt
        sub!(rop, x, y)
    elseif _fits_clong(y)
        sub!(rop, x, _as_clong(y))
    elseif _fits_culong(y)
        sub!(rop, x, _as_culong(y))
    else
        sub!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function sub!(rop::BigInt, x::Integer, y::BigInt)
    if x isa NativeInt
        sub!(rop, x, y)
    elseif _fits_clong(x)
        sub!(rop, _as_clong(x), y)
    elseif _fits_culong(x)
        sub!(rop, _as_culong(x), y)
    else
        sub!(rop, _as_bigint(x), y)
    end
    return rop
end

@inline function sub!(rop::BigInt, x::T, y::U) where {T <: NativeInt, U <: NativeInt}
    return set!(rop, BigInt(x) - BigInt(y))
end

@inline function sub!(rop::BigInt, x::Integer, y::Integer)
    return set!(rop, _as_bigint(x) - _as_bigint(y))
end

# ============================================================================
# BigInt: mul! (Unified multi-dispatch with mpz_mul_ui optimization)
# ============================================================================

@inline function mul!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_mul(rop, x, y)
    return rop
end

@inline function mul!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_mul_ui(rop, x, y)
    return rop
end

@inline function mul!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    if y < 0
        mpz_mul_ui(rop, x, Culong(-y))
        mpz_neg(rop, rop)
    else
        mpz_mul_ui(rop, x, Culong(y))
    end
    return rop
end

@inline function mul!(rop::BigInt, x::T, y::BigInt) where {T <: NativeInt}
    return mul!(rop, y, x)
end

@inline function mul!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeInt
        mul!(rop, x, y)
    elseif _fits_clong(y)
        mul!(rop, x, _as_clong(y))
    elseif _fits_culong(y)
        mul!(rop, x, _as_culong(y))
    else
        mul!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function mul!(rop::BigInt, x::Integer, y::BigInt)
    if x isa NativeInt
        mul!(rop, y, x)
    elseif _fits_clong(x)
        mul!(rop, y, _as_clong(x))
    elseif _fits_culong(x)
        mul!(rop, y, _as_culong(x))
    else
        mul!(rop, _as_bigint(x), y)
    end
    return rop
end

@inline function mul!(rop::BigInt, x::T, y::U) where {T <: NativeInt, U <: NativeInt}
    return set!(rop, BigInt(x) * BigInt(y))
end

@inline function mul!(rop::BigInt, x::Integer, y::Integer)
    return set!(rop, _as_bigint(x) * _as_bigint(y))
end

# ============================================================================
# BigInt: Compound operations (addmul!, submul!, mul_2exp!)
# ============================================================================

@inline function addmul!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_addmul(rop, x, y)
    return rop
end

@inline function addmul!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_addmul_ui(rop, x, y)
    return rop
end

@inline function addmul!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    if y < 0
        mpz_submul_ui(rop, x, Culong(-y))
    else
        mpz_addmul_ui(rop, x, Culong(y))
    end
    return rop
end

@inline function addmul!(rop::BigInt, x::Integer, y::BigInt)
    return addmul!(rop, y, x)
end

@inline function addmul!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeInt
        addmul!(rop, x, y)
    elseif _fits_clong(y)
        addmul!(rop, x, _as_clong(y))
    elseif _fits_culong(y)
        addmul!(rop, x, _as_culong(y))
    else
        addmul!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function addmul!(rop::BigInt, x::Integer, y::Integer)
    addmul!(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

@inline function submul!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_submul(rop, x, y)
    return rop
end

@inline function submul!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_submul_ui(rop, x, y)
    return rop
end

@inline function submul!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    if y < 0
        mpz_addmul_ui(rop, x, Culong(-y))
    else
        mpz_submul_ui(rop, x, Culong(y))
    end
    return rop
end

@inline function submul!(rop::BigInt, x::Integer, y::BigInt)
    return submul!(rop, y, x)
end

@inline function submul!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeInt
        submul!(rop, x, y)
    elseif _fits_clong(y)
        submul!(rop, x, _as_clong(y))
    elseif _fits_culong(y)
        submul!(rop, x, _as_culong(y))
    else
        submul!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function submul!(rop::BigInt, x::Integer, y::Integer)
    submul!(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

@inline function mul_2exp!(rop::BigInt, x::BigInt, n::T) where {T <: NativeInt}
    _checked_nonnegative(Val(:mul_2exp!), n)
    mpz_mul_2exp(rop, x, n)
    return rop
end

@inline function mul_2exp!(rop::BigInt, x::Integer, n::Integer)
    _checked_nonnegative(Val(:mul_2exp!), n)
    mpz_mul_2exp(rop, _as_bigint(x), n)
    return rop
end

@inline function div_2exp!(rop::BigInt, x::BigInt, n::T; rounding::Symbol=:trunc) where {T <: NativeInt}
    _checked_nonnegative(Val(:div_2exp!), n)
    if rounding === :trunc
        mpz_tdiv_q_2exp(rop, x, n)
    elseif rounding === :floor
        mpz_fdiv_q_2exp(rop, x, n)
    elseif rounding === :ceil
        mpz_cdiv_q_2exp(rop, x, n)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return rop
end

@inline function div_2exp!(rop::BigInt, x::Integer, n::Integer; rounding::Symbol=:trunc)
    _checked_nonnegative(Val(:div_2exp!), n)
    if rounding === :trunc
        mpz_tdiv_q_2exp(rop, _as_bigint(x), n)
    elseif rounding === :floor
        mpz_fdiv_q_2exp(rop, _as_bigint(x), n)
    elseif rounding === :ceil
        mpz_cdiv_q_2exp(rop, _as_bigint(x), n)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return rop
end

# ============================================================================
# BigInt: Bit manipulation (setbit!, clrbit!, combit!)
# ============================================================================

@inline function setbit!(rop::BigInt, bit::T) where {T <: NativeInt}
    _checked_nonnegative(Val(:setbit!), bit)
    mpz_setbit(rop, bit)
    return rop
end

@inline function setbit!(rop::BigInt, bit::Integer)
    _checked_nonnegative(Val(:setbit!), bit)
    mpz_setbit(rop, bit)
    return rop
end

@inline function clrbit!(rop::BigInt, bit::T) where {T <: NativeInt}
    _checked_nonnegative(Val(:clrbit!), bit)
    mpz_clrbit(rop, bit)
    return rop
end

@inline function clrbit!(rop::BigInt, bit::Integer)
    _checked_nonnegative(Val(:clrbit!), bit)
    mpz_clrbit(rop, bit)
    return rop
end

@inline function combit!(rop::BigInt, bit::T) where {T <: NativeInt}
    _checked_nonnegative(Val(:combit!), bit)
    mpz_combit(rop, bit)
    return rop
end

@inline function combit!(rop::BigInt, bit::Integer)
    _checked_nonnegative(Val(:combit!), bit)
    mpz_combit(rop, bit)
    return rop
end

# ============================================================================
# BigInt: Division operations (unified with dispatch)
# ============================================================================

@inline function div!(rop::BigInt, x::BigInt, y::BigInt; rounding::Symbol=:trunc)
    if rounding === :trunc
        mpz_tdiv_q(rop, x, y)
    elseif rounding === :floor
        mpz_fdiv_q(rop, x, y)
    elseif rounding === :ceil
        mpz_cdiv_q(rop, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return rop
end

@inline function div!(rop::BigInt, x::BigInt, y::T; rounding::Symbol=:trunc) where {T <: NativeUnsignedInt}
    if rounding === :trunc
        mpz_tdiv_q_ui(rop, x, y)
    elseif rounding === :floor
        mpz_fdiv_q_ui(rop, x, y)
    elseif rounding === :ceil
        mpz_cdiv_q_ui(rop, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return rop
end

@inline function div!(rop::BigInt, x::BigInt, y::T; rounding::Symbol=:trunc) where {T <: NativeSignedInt}
    return div!(rop, x, _as_bigint(y); rounding=rounding)
end

@inline function div!(rop::BigInt, x::Integer, y::BigInt; rounding::Symbol=:trunc)
    return div!(rop, _as_bigint(x), y; rounding=rounding)
end

@inline function div!(rop::BigInt, x::BigInt, y::Integer; rounding::Symbol=:trunc)
    if y isa NativeUnsignedInt
        div!(rop, x, y; rounding=rounding)
    elseif _fits_culong(y)
        div!(rop, x, _as_culong(y); rounding=rounding)
    else
        div!(rop, x, _as_bigint(y); rounding=rounding)
    end
    return rop
end

@inline function div!(rop::BigInt, x::Integer, y::Integer; rounding::Symbol=:trunc)
    return div!(rop, _as_bigint(x), _as_bigint(y); rounding=rounding)
end

@inline function rem!(rop::BigInt, x::BigInt, y::BigInt; rounding::Symbol=:trunc)
    if rounding === :trunc
        mpz_tdiv_r(rop, x, y)
    elseif rounding === :floor
        mpz_fdiv_r(rop, x, y)
    elseif rounding === :ceil
        mpz_cdiv_r(rop, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return rop
end

@inline function rem!(rop::BigInt, x::BigInt, y::T; rounding::Symbol=:trunc) where {T <: NativeUnsignedInt}
    if rounding === :trunc
        mpz_tdiv_r_ui(rop, x, y)
    elseif rounding === :floor
        mpz_fdiv_r_ui(rop, x, y)
    elseif rounding === :ceil
        mpz_cdiv_r_ui(rop, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return rop
end

@inline function rem!(rop::BigInt, x::BigInt, y::T; rounding::Symbol=:trunc) where {T <: NativeSignedInt}
    return rem!(rop, x, _as_bigint(y); rounding=rounding)
end

@inline function rem!(rop::BigInt, x::Integer, y::BigInt; rounding::Symbol=:trunc)
    return rem!(rop, _as_bigint(x), y; rounding=rounding)
end

@inline function rem!(rop::BigInt, x::BigInt, y::Integer; rounding::Symbol=:trunc)
    if y isa NativeUnsignedInt
        rem!(rop, x, y; rounding=rounding)
    elseif _fits_culong(y)
        rem!(rop, x, _as_culong(y); rounding=rounding)
    else
        rem!(rop, x, _as_bigint(y); rounding=rounding)
    end
    return rop
end

@inline function rem!(rop::BigInt, x::Integer, y::Integer; rounding::Symbol=:trunc)
    return rem!(rop, _as_bigint(x), _as_bigint(y); rounding=rounding)
end

@inline function divrem!(q::BigInt, r::BigInt, x::BigInt, y::BigInt; rounding::Symbol=:trunc)
    if rounding === :trunc
        mpz_tdiv_qr(q, r, x, y)
    elseif rounding === :floor
        mpz_fdiv_qr(q, r, x, y)
    elseif rounding === :ceil
        mpz_cdiv_qr(q, r, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return q
end

@inline function divrem!(q::BigInt, r::BigInt, x::BigInt, y::T; rounding::Symbol=:trunc) where {T <: NativeUnsignedInt}
    if rounding === :trunc
        mpz_tdiv_qr_ui(q, r, x, y)
    elseif rounding === :floor
        mpz_fdiv_qr_ui(q, r, x, y)
    elseif rounding === :ceil
        mpz_cdiv_qr_ui(q, r, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return q
end

@inline function divrem!(q::BigInt, r::BigInt, x::Integer, y::Integer; rounding::Symbol=:trunc)
    return divrem!(q, r, _as_bigint(x), _as_bigint(y); rounding=rounding)
end

@inline function divrem_ui!(q::BigInt, r::Ref{Culong}, x::BigInt, y::T; rounding::Symbol=:trunc) where {T <: NativeUnsignedInt}
    if rounding === :trunc
        mpz_tdiv_qr_ui(q, r, x, y)
    elseif rounding === :floor
        mpz_fdiv_qr_ui(q, r, x, y)
    elseif rounding === :ceil
        mpz_cdiv_qr_ui(q, r, x, y)
    else
        throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    end
    return q
end

@inline function mod!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_mod(rop, x, y)
    return rop
end

@inline function mod!(rop::BigInt, x::Integer, y::BigInt)
    return mod!(rop, _as_bigint(x), y)
end

@inline function mod!(rop::BigInt, x::BigInt, y::Integer)
    return mod!(rop, x, _as_bigint(y))
end

@inline function mod!(rop::BigInt, x::Integer, y::Integer)
    return mod!(rop, _as_bigint(x), _as_bigint(y))
end

@inline function divexact!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_divexact(rop, x, y)
    return rop
end

@inline function divexact!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_divexact_ui(rop, x, y)
    return rop
end

@inline function divexact!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    return divexact!(rop, x, _as_bigint(y))
end

@inline function divexact!(rop::BigInt, x::Integer, y::BigInt)
    return divexact!(rop, _as_bigint(x), y)
end

@inline function divexact!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeUnsignedInt
        divexact!(rop, x, y)
    elseif _fits_culong(y)
        divexact!(rop, x, _as_culong(y))
    else
        divexact!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function divexact!(rop::BigInt, x::Integer, y::Integer)
    return divexact!(rop, _as_bigint(x), _as_bigint(y))
end

# ============================================================================
# BigInt: Power functions (pow!, powm!, etc.)
# ============================================================================

@inline function pow!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    _checked_nonnegative(Val(:pow!), y)
    mpz_pow_ui(rop, x, y)
    return rop
end

@inline function pow!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    _checked_nonnegative(Val(:pow!), y)
    mpz_pow_ui(rop, x, Culong(y))
    return rop
end

@inline function pow!(rop::BigInt, x::T, y::T) where {T <: NativeUnsignedInt}
    _checked_nonnegative(Val(:pow!), y)
    mpz_ui_pow_ui(rop, x, y)
    return rop
end

@inline function pow!(rop::BigInt, x::T, y::U) where {T <: NativeSignedInt, U <: NativeInt}
    _checked_nonnegative(Val(:pow!), y)
    if x < 0
        mpz_ui_pow_ui(rop, Culong(-x), Culong(y))
        isodd(y) && mpz_neg(rop, rop)
    else
        mpz_ui_pow_ui(rop, Culong(x), Culong(y))
    end
    return rop
end

@inline function pow!(rop::BigInt, x::Integer, y::Integer)
    _checked_nonnegative(Val(:pow!), y)
    if y <= typemax(Culong)
        xb = _as_bigint(x)
        mpz_pow_ui(rop, xb, Culong(y))
    else
        set!(rop, _as_bigint(x) ^ y)
    end
    return rop
end

@inline function powm!(rop::BigInt, x::BigInt, y::BigInt, m::BigInt)
    mpz_powm(rop, x, y, m)
    return rop
end

@inline function powm!(rop::BigInt, x::BigInt, y::T, m::BigInt) where {T <: NativeUnsignedInt}
    mpz_powm_ui(rop, x, y, m)
    return rop
end

@inline function powm!(rop::BigInt, x::Integer, y::Integer, m::BigInt)
    return powm!(rop, _as_bigint(x), _as_bigint(y), m)
end

@inline function powm_sec!(rop::BigInt, x::BigInt, y::BigInt, m::BigInt)
    mpz_powm_sec(rop, x, y, m)
    return rop
end

@inline function powm_sec!(rop::BigInt, x::Integer, y::Integer, m::BigInt)
    return powm_sec!(rop, _as_bigint(x), _as_bigint(y), m)
end

# ============================================================================
# BigInt: Root functions (sqrt!, root!, rootrem!, sqrtrem!)
# ============================================================================

@inline function sqrt!(rop::BigInt, x::BigInt)
    mpz_sqrt(rop, x)
    return rop
end

@inline function sqrt!(rop::BigInt, x::Integer)
    return sqrt!(rop, _as_bigint(x))
end

@inline function sqrtrem!(rop1::BigInt, rop2::BigInt, x::BigInt)
    mpz_sqrtrem(rop1, rop2, x)
    return rop1
end

@inline function sqrtrem!(rop1::BigInt, rop2::BigInt, x::Integer)
    return sqrtrem!(rop1, rop2, _as_bigint(x))
end

@inline function root!(rop::BigInt, x::BigInt, n::T) where {T <: NativeUnsignedInt}
    return mpz_root(rop, x, n)
end

@inline function root!(rop::BigInt, x::Integer, n::Integer)
    _checked_nonnegative(Val(:root!), n)
    n <= typemax(Culong) || throw(ArgumentError("root! exponent does not fit in Culong"))
    return mpz_root(rop, _as_bigint(x), Culong(n))
end

@inline function rootrem!(rop1::BigInt, rop2::BigInt, x::BigInt, n::T) where {T <: NativeUnsignedInt}
    mpz_rootrem(rop1, rop2, x, n)
    return rop1
end

@inline function rootrem!(rop1::BigInt, rop2::BigInt, x::Integer, n::Integer)
    _checked_nonnegative(Val(:rootrem!), n)
    n <= typemax(Culong) || throw(ArgumentError("rootrem! exponent does not fit in Culong"))
    mpz_rootrem(rop1, rop2, _as_bigint(x), Culong(n))
    return rop1
end

# ============================================================================
# BigInt: Unary operations (neg!, abs!)
# ============================================================================

@inline function neg!(rop::BigInt, x::BigInt)
    mpz_neg(rop, x)
    return rop
end

@inline function neg!(rop::BigInt, x::Integer)
    return neg!(rop, _as_bigint(x))
end

@inline function abs!(rop::BigInt, x::BigInt)
    mpz_abs(rop, x)
    return rop
end

@inline function abs!(rop::BigInt, x::Integer)
    return abs!(rop, _as_bigint(x))
end

# ============================================================================
# BigInt: Number Theoretic Functions (gcd, lcm, invert, remove, etc.)
# ============================================================================

@inline function gcd!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_gcd(rop, x, y)
    return rop
end

@inline function gcd!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    return mpz_gcd_ui(rop, x, y)
end

@inline function gcd!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    return gcd!(rop, x, _as_bigint(y))
end

@inline function gcd!(rop::BigInt, x::Integer, y::BigInt)
    return gcd!(rop, _as_bigint(x), y)
end

@inline function gcd!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeUnsignedInt
        gcd!(rop, x, y)
    else
        gcd!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function gcd!(rop::BigInt, x::Integer, y::Integer)
    return gcd!(rop, _as_bigint(x), _as_bigint(y))
end

@inline function gcd_ui!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    return mpz_gcd_ui(rop, x, y)
end

@inline function gcdext!(rop::BigInt, s::BigInt, t::BigInt, x::BigInt, y::BigInt)
    mpz_gcdext(rop, s, t, x, y)
    return rop
end

@inline function gcdext!(rop::BigInt, s::BigInt, t::BigInt, x::Integer, y::Integer)
    return gcdext!(rop, s, t, _as_bigint(x), _as_bigint(y))
end

@inline function lcm!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_lcm(rop, x, y)
    return rop
end

@inline function lcm!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_lcm_ui(rop, x, y)
    return rop
end

@inline function lcm!(rop::BigInt, x::BigInt, y::T) where {T <: NativeSignedInt}
    return lcm!(rop, x, _as_bigint(y))
end

@inline function lcm!(rop::BigInt, x::Integer, y::BigInt)
    return lcm!(rop, _as_bigint(x), y)
end

@inline function lcm!(rop::BigInt, x::BigInt, y::Integer)
    if y isa NativeUnsignedInt
        lcm!(rop, x, y)
    else
        lcm!(rop, x, _as_bigint(y))
    end
    return rop
end

@inline function lcm!(rop::BigInt, x::Integer, y::Integer)
    return lcm!(rop, _as_bigint(x), _as_bigint(y))
end

@inline function lcm_ui!(rop::BigInt, x::BigInt, y::T) where {T <: NativeUnsignedInt}
    mpz_lcm_ui(rop, x, y)
    return rop
end

@inline function invert!(rop::BigInt, x::BigInt, m::BigInt)
    rc = mpz_invert(rop, x, m)
    rc != 0 || error("inverse does not exist")
    return rop
end

@inline function invert!(rop::BigInt, x::Integer, m::BigInt)
    return invert!(rop, _as_bigint(x), m)
end

@inline function remove!(rop::BigInt, x::BigInt, f::BigInt)
    return mpz_remove(rop, x, f)
end

@inline function remove!(rop::BigInt, x::Integer, f::Integer)
    return remove!(rop, _as_bigint(x), _as_bigint(f))
end

# ============================================================================
# BigInt: Combinatorial functions (fac_ui!, bin_ui!, fib_ui!, lucnum_ui!, etc.)
# ============================================================================

@inline function fac_ui!(rop::BigInt, n::T) where {T <: NativeUnsignedInt}
    mpz_fac_ui(rop, n)
    return rop
end

@inline function fac_ui!(rop::BigInt, n::Integer)
    _checked_nonnegative(Val(:fac_ui!), n)
    n <= typemax(Culong) || throw(ArgumentError("fac_ui! argument does not fit in Culong"))
    mpz_fac_ui(rop, Culong(n))
    return rop
end

@inline function bin_ui!(rop::BigInt, n::BigInt, k::T) where {T <: NativeUnsignedInt}
    mpz_bin_ui(rop, n, k)
    return rop
end

@inline function bin_ui!(rop::BigInt, n::Integer, k::Integer)
    _checked_nonnegative(Val(:bin_ui!), k)
    k <= typemax(Culong) || throw(ArgumentError("bin_ui! argument does not fit in Culong"))
    mpz_bin_ui(rop, _as_bigint(n), Culong(k))
    return rop
end

@inline function bin_uiui!(rop::BigInt, n::T, k::T) where {T <: NativeUnsignedInt}
    mpz_bin_uiui(rop, n, k)
    return rop
end

@inline function bin_uiui!(rop::BigInt, n::Integer, k::Integer)
    _checked_nonnegative(Val(:bin_uiui!), n)
    _checked_nonnegative(Val(:bin_uiui!), k)
    n <= typemax(Culong) || throw(ArgumentError("bin_uiui! argument does not fit in Culong"))
    k <= typemax(Culong) || throw(ArgumentError("bin_uiui! argument does not fit in Culong"))
    mpz_bin_uiui(rop, Culong(n), Culong(k))
    return rop
end

@inline function fib_ui!(rop::BigInt, n::T) where {T <: NativeUnsignedInt}
    mpz_fib_ui(rop, n)
    return rop
end

@inline function fib_ui!(rop::BigInt, n::Integer)
    _checked_nonnegative(Val(:fib_ui!), n)
    n <= typemax(Culong) || throw(ArgumentError("fib_ui! argument does not fit in Culong"))
    mpz_fib_ui(rop, Culong(n))
    return rop
end

@inline function fib2_ui!(rop1::BigInt, rop2::BigInt, n::T) where {T <: NativeUnsignedInt}
    mpz_fib2_ui(rop1, rop2, n)
    return rop1
end

@inline function fib2_ui!(rop1::BigInt, rop2::BigInt, n::Integer)
    _checked_nonnegative(Val(:fib2_ui!), n)
    n <= typemax(Culong) || throw(ArgumentError("fib2_ui! argument does not fit in Culong"))
    mpz_fib2_ui(rop1, rop2, Culong(n))
    return rop1
end

@inline function lucnum_ui!(rop::BigInt, n::T) where {T <: NativeUnsignedInt}
    mpz_lucnum_ui(rop, n)
    return rop
end

@inline function lucnum_ui!(rop::BigInt, n::Integer)
    _checked_nonnegative(Val(:lucnum_ui!), n)
    n <= typemax(Culong) || throw(ArgumentError("lucnum_ui! argument does not fit in Culong"))
    mpz_lucnum_ui(rop, Culong(n))
    return rop
end

@inline function lucnum2_ui!(rop1::BigInt, rop2::BigInt, n::T) where {T <: NativeUnsignedInt}
    mpz_lucnum2_ui(rop1, rop2, n)
    return rop1
end

@inline function lucnum2_ui!(rop1::BigInt, rop2::BigInt, n::Integer)
    _checked_nonnegative(Val(:lucnum2_ui!), n)
    n <= typemax(Culong) || throw(ArgumentError("lucnum2_ui! argument does not fit in Culong"))
    mpz_lucnum2_ui(rop1, rop2, Culong(n))
    return rop1
end

# ============================================================================
# BigInt: Bitwise operations (and!, or!, xor!, com!)
# ============================================================================

@inline function and!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_and(rop, x, y)
    return rop
end

@inline function and!(rop::BigInt, x::Integer, y::BigInt)
    return and!(rop, _as_bigint(x), y)
end

@inline function and!(rop::BigInt, x::BigInt, y::Integer)
    return and!(rop, x, _as_bigint(y))
end

@inline function and!(rop::BigInt, x::Integer, y::Integer)
    return and!(rop, _as_bigint(x), _as_bigint(y))
end

@inline function or!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_ior(rop, x, y)
    return rop
end

@inline function or!(rop::BigInt, x::Integer, y::BigInt)
    return or!(rop, _as_bigint(x), y)
end

@inline function or!(rop::BigInt, x::BigInt, y::Integer)
    return or!(rop, x, _as_bigint(y))
end

@inline function or!(rop::BigInt, x::Integer, y::Integer)
    return or!(rop, _as_bigint(x), _as_bigint(y))
end

@inline function xor!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_xor(rop, x, y)
    return rop
end

@inline function xor!(rop::BigInt, x::Integer, y::BigInt)
    return xor!(rop, _as_bigint(x), y)
end

@inline function xor!(rop::BigInt, x::BigInt, y::Integer)
    return xor!(rop, x, _as_bigint(y))
end

@inline function xor!(rop::BigInt, x::Integer, y::Integer)
    return xor!(rop, _as_bigint(x), _as_bigint(y))
end

@inline function com!(rop::BigInt, x::BigInt)
    mpz_com(rop, x)
    return rop
end

@inline function com!(rop::BigInt, x::Integer)
    return com!(rop, _as_bigint(x))
end

@inline function nextprime!(rop::BigInt, x::BigInt)
    mpz_nextprime(rop, x)
    return rop
end

@inline function nextprime!(rop::BigInt, x::Integer)
    return nextprime!(rop, _as_bigint(x))
end

# ============================================================================
# BigFloat: Assignment and conversion (Tier 1/2 multi-dispatch)
# ============================================================================

@inline function set!(rop::BigFloat, x::BigFloat; rounding=MPFRRoundNearest)
    mpfr_set(rop, x, _rounding(rounding))
    return rop
end

@inline function set!(rop::BigFloat, x::BigInt; rounding=MPFRRoundNearest)
    mpfr_set_z(rop, x, _rounding(rounding))
    return rop
end

@inline function set!(rop::BigFloat, x::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_set_si(rop, x, _rounding(rounding))
    return rop
end

@inline function set!(rop::BigFloat, x::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_set_ui(rop, x, _rounding(rounding))
    return rop
end

@inline function set!(rop::BigFloat, x::AbstractFloat; rounding=MPFRRoundNearest)
    mpfr_set_d(rop, x, _rounding(rounding))
    return rop
end

@inline function set!(rop::BigFloat, x::Integer; rounding=MPFRRoundNearest)
    if x isa BigInt
        mpfr_set_z(rop, x, _rounding(rounding))
    elseif x isa NativeSignedInt
        mpfr_set_si(rop, x, _rounding(rounding))
    elseif x isa NativeUnsignedInt
        mpfr_set_ui(rop, x, _rounding(rounding))
    else
        set!(rop, BigInt(x); rounding=rounding)
    end
    return rop
end

function set!(rop::BigFloat, s::AbstractString; base::Integer=0, rounding=MPFRRoundNearest)
    rc = mpfr_set_str(rop, s, base, _rounding(rounding))
    rc == 0 || throw(ArgumentError("invalid floating-point string for base $base: $s"))
    return rop
end

@inline function swap!(x::BigFloat, y::BigFloat)
    tmp = BigFloat(0; precision=precision(x))
    set!(tmp, x)
    set!(x, y)
    set!(y, tmp)
    return x
end

# ============================================================================
# BigFloat: Binary Arithmetic (add!, sub!, mul!, div!) - Unified dispatch
# ============================================================================

@inline function add!(rop::BigFloat, x::BigFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_add(rop, x, y, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::BigFloat, y::BigInt; rounding=MPFRRoundNearest)
    mpfr_add_z(rop, x, y, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::BigInt, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_add_z(rop, y, x, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_add_ui(rop, x, y, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_add_ui(rop, y, x, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_add_si(rop, x, y, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_add_si(rop, y, x, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::BigFloat, y::AbstractFloat; rounding=MPFRRoundNearest)
    mpfr_add_d(rop, x, y, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::AbstractFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_add_d(rop, y, x, _rounding(rounding))
    return rop
end

@inline function add!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    if x isa BigInt && y isa BigInt
        mpfr_add(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    elseif x isa BigInt
        add!(rop, y, x; rounding=r)
    elseif y isa BigInt
        add!(rop, x, y; rounding=r)
    else
        mpfr_add(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    end
    return rop
end

@inline function sub!(rop::BigFloat, x::BigFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_sub(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::BigFloat, y::BigInt; rounding=MPFRRoundNearest)
    mpfr_sub_z(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::BigInt, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_z_sub(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_sub_ui(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_ui_sub(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_sub_si(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_si_sub(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::BigFloat, y::AbstractFloat; rounding=MPFRRoundNearest)
    mpfr_sub_d(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::AbstractFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_d_sub(rop, x, y, _rounding(rounding))
    return rop
end

@inline function sub!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    if x isa BigInt && y isa BigInt
        mpfr_sub(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    elseif x isa BigInt
        mpfr_z_sub(rop, x, _as_bigfloat(y; precision=precision(rop)), r)
    elseif y isa BigInt
        mpfr_sub_z(rop, _as_bigfloat(x; precision=precision(rop)), y, r)
    else
        mpfr_sub(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    end
    return rop
end

@inline function mul!(rop::BigFloat, x::BigFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_mul(rop, x, y, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::BigFloat, y::BigInt; rounding=MPFRRoundNearest)
    mpfr_mul_z(rop, x, y, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::BigInt, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_mul_z(rop, y, x, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_mul_ui(rop, x, y, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_mul_ui(rop, y, x, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_mul_si(rop, x, y, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_mul_si(rop, y, x, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::BigFloat, y::AbstractFloat; rounding=MPFRRoundNearest)
    mpfr_mul_d(rop, x, y, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::AbstractFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_mul_d(rop, y, x, _rounding(rounding))
    return rop
end

@inline function mul!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    if x isa BigInt && y isa BigInt
        mpfr_mul(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    elseif x isa BigInt
        mpfr_mul_z(rop, _as_bigfloat(y; precision=precision(rop)), x, r)
    elseif y isa BigInt
        mpfr_mul_z(rop, _as_bigfloat(x; precision=precision(rop)), y, r)
    else
        mpfr_mul(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    end
    return rop
end

@inline function div!(rop::BigFloat, x::BigFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_div(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::BigFloat, y::BigInt; rounding=MPFRRoundNearest)
    mpfr_div_z(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::BigInt, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_div(rop, _as_bigfloat(x; precision=precision(rop)), y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_div_ui(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_ui_div(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::BigFloat, y::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_div_si(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::T, y::BigFloat; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_si_div(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::BigFloat, y::AbstractFloat; rounding=MPFRRoundNearest)
    mpfr_div_d(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::AbstractFloat, y::BigFloat; rounding=MPFRRoundNearest)
    mpfr_d_div(rop, x, y, _rounding(rounding))
    return rop
end

@inline function div!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    if x isa BigInt && y isa BigInt
        mpfr_div(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    elseif x isa BigInt
        mpfr_div(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    elseif y isa BigInt
        mpfr_div_z(rop, _as_bigfloat(x; precision=precision(rop)), y, r)
    else
        mpfr_div(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    end
    return rop
end

# ============================================================================
# BigFloat: Unary operations (neg!, abs!, sqrt!, square!)
# ============================================================================

@inline function neg!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_neg(z, a, r), rop, x; rounding=rounding)
end

@inline function abs!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_abs(z, a, r), rop, x; rounding=rounding)
end

@inline function sqrt!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sqrt(z, a, r), rop, x; rounding=rounding)
end

@inline function square!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sqr(z, a, r), rop, x; rounding=rounding)
end

@inline function rec_sqrt!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_rec_sqrt(z, a, r), rop, x; rounding=rounding)
end

@inline function cbrt!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_cbrt(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Root extraction - Unified dispatch (rootn! handles ui/si internally)
# ============================================================================

@inline function rootn!(rop::BigFloat, x::FloatLike, n::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    p = precision(rop)
    mpfr_rootn_ui(rop, _as_bigfloat(x; precision=p), n, _rounding(rounding))
    return rop
end

@inline function rootn!(rop::BigFloat, x::FloatLike, n::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    p = precision(rop)
    mpfr_rootn_si(rop, _as_bigfloat(x; precision=p), n, _rounding(rounding))
    return rop
end

# ============================================================================
# BigFloat: Power functions - Unified dispatch (pow!, pow_z)
# ============================================================================

@inline function pow!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_pow(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function pow!(rop::BigFloat, x::FloatLike, y::BigInt; rounding=MPFRRoundNearest)
    mpfr_pow_z(rop, _as_bigfloat(x; precision=precision(rop)), y, _rounding(rounding))
    return rop
end

@inline function pow!(rop::BigFloat, x::FloatLike, y::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_pow_ui(rop, _as_bigfloat(x; precision=precision(rop)), y, _rounding(rounding))
    return rop
end

@inline function pow!(rop::BigFloat, x::FloatLike, y::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_pow_si(rop, _as_bigfloat(x; precision=precision(rop)), y, _rounding(rounding))
    return rop
end

@inline function pow!(rop::BigFloat, x::T, y::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_ui_pow_ui(rop, x, y, _rounding(rounding))
    return rop
end

# Alias for consistency
@inline function pow!(rop::BigFloat, x::T, y::FloatLike; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_ui_pow(rop, x, _as_bigfloat(y; precision=precision(rop)), _rounding(rounding))
    return rop
end

# ============================================================================
# BigFloat: Exponential functions (exp, exp2, exp10, exp1m, etc.) - Unified
# ============================================================================

@inline function exp!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_exp(z, a, r), rop, x; rounding=rounding)
end

@inline function exp2!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_exp2(z, a, r), rop, x; rounding=rounding)
end

@inline function exp10!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_exp10(z, a, r), rop, x; rounding=rounding)
end

@inline function expm1!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_expm1(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Logarithmic functions (log, log2, log10, log1p) - Unified
# ============================================================================

@inline function log!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_log(z, a, r), rop, x; rounding=rounding)
end

@inline function log2!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_log2(z, a, r), rop, x; rounding=rounding)
end

@inline function log10!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_log10(z, a, r), rop, x; rounding=rounding)
end

@inline function log1p!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_log1p(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Trigonometric functions (sin, cos, tan, sin_cos, etc.) - Unified
# ============================================================================

@inline function sin!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sin(z, a, r), rop, x; rounding=rounding)
end

@inline function cos!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_cos(z, a, r), rop, x; rounding=rounding)
end

@inline function tan!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_tan(z, a, r), rop, x; rounding=rounding)
end

@inline function sin_cos!(rop_sin::BigFloat, rop_cos::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    p = precision(rop_sin)
    mpfr_sin_cos(rop_sin, rop_cos, _as_bigfloat(x; precision=p), _rounding(rounding))
    return rop_sin
end

@inline function sec!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sec(z, a, r), rop, x; rounding=rounding)
end

@inline function csc!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_csc(z, a, r), rop, x; rounding=rounding)
end

@inline function cot!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_cot(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Inverse trigonometric functions (acos, asin, atan, atan2) - Unified
# ============================================================================

@inline function acos!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_acos(z, a, r), rop, x; rounding=rounding)
end

@inline function asin!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_asin(z, a, r), rop, x; rounding=rounding)
end

@inline function atan!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_atan(z, a, r), rop, x; rounding=rounding)
end

@inline function atan2!(rop::BigFloat, y::FloatLike, x::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_atan2(rop, _as_bigfloat(y; precision=p), _as_bigfloat(x; precision=p), r)
    return rop
end

# ============================================================================
# BigFloat: Hyperbolic functions (sinh, cosh, tanh, sinh_cosh, etc.) - Unified
# ============================================================================

@inline function sinh!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sinh(z, a, r), rop, x; rounding=rounding)
end

@inline function cosh!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_cosh(z, a, r), rop, x; rounding=rounding)
end

@inline function tanh!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_tanh(z, a, r), rop, x; rounding=rounding)
end

@inline function sinh_cosh!(rop_sinh::BigFloat, rop_cosh::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    p = precision(rop_sinh)
    mpfr_sinh_cosh(rop_sinh, rop_cosh, _as_bigfloat(x; precision=p), _rounding(rounding))
    return rop_sinh
end

@inline function sech!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sech(z, a, r), rop, x; rounding=rounding)
end

@inline function csch!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_csch(z, a, r), rop, x; rounding=rounding)
end

@inline function coth!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_coth(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Inverse hyperbolic functions (acosh, asinh, atanh) - Unified
# ============================================================================

@inline function acosh!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_acosh(z, a, r), rop, x; rounding=rounding)
end

@inline function asinh!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_asinh(z, a, r), rop, x; rounding=rounding)
end

@inline function atanh!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_atanh(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Special functions (gamma, lgamma, digamma, beta, zeta, erf, erfc, etc.)
# ============================================================================

@inline function gamma!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_gamma(z, a, r), rop, x; rounding=rounding)
end

@inline function gamma_inc!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_gamma_inc(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function lngamma!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_lngamma(z, a, r), rop, x; rounding=rounding)
end

@inline function lgamma!(rop::BigFloat, signp::Ref{Cint}, x::FloatLike; rounding=MPFRRoundNearest)
    p = precision(rop)
    mpfr_lgamma(rop, signp, _as_bigfloat(x; precision=p), _rounding(rounding))
    return rop
end

@inline function digamma!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_digamma(z, a, r), rop, x; rounding=rounding)
end

@inline function beta!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_beta(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function zeta!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_zeta(z, a, r), rop, x; rounding=rounding)
end

@inline function zeta!(rop::BigFloat, n::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_zeta_ui(rop, n, _rounding(rounding))
    return rop
end

@inline function erf!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_erf(z, a, r), rop, x; rounding=rounding)
end

@inline function erfc!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_erfc(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Bessel functions (j0, j1, jn, y0, y1, yn) - Unified with dispatch
# ============================================================================

@inline function j0!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_j0(z, a, r), rop, x; rounding=rounding)
end

@inline function j1!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_j1(z, a, r), rop, x; rounding=rounding)
end

@inline function jn!(rop::BigFloat, n::T, x::FloatLike; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_jn(rop, n, _as_bigfloat(x; precision=precision(rop)), _rounding(rounding))
    return rop
end

@inline function y0!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_y0(z, a, r), rop, x; rounding=rounding)
end

@inline function y1!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_y1(z, a, r), rop, x; rounding=rounding)
end

@inline function yn!(rop::BigFloat, n::T, x::FloatLike; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_yn(rop, n, _as_bigfloat(x; precision=precision(rop)), _rounding(rounding))
    return rop
end

# ============================================================================
# BigFloat: Other special functions (eint, li2, agm, ai)
# ============================================================================

@inline function eint!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_eint(z, a, r), rop, x; rounding=rounding)
end

@inline function li2!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_li2(z, a, r), rop, x; rounding=rounding)
end

@inline function agm!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_agm(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function ai!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_ai(z, a, r), rop, x; rounding=rounding)
end

# ============================================================================
# BigFloat: Rounding and remainder functions (rint, ceil, floor, round, trunc, frac, modf)
# ============================================================================

@inline function rint!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_rint(z, a, r), rop, x; rounding=rounding)
end

@inline function ceil!(rop::BigFloat, x::FloatLike)
    p = precision(rop)
    mpfr_ceil(rop, _as_bigfloat(x; precision=p))
    return rop
end

@inline function floor!(rop::BigFloat, x::FloatLike)
    p = precision(rop)
    mpfr_floor(rop, _as_bigfloat(x; precision=p))
    return rop
end

@inline function round!(rop::BigFloat, x::FloatLike)
    p = precision(rop)
    mpfr_round(rop, _as_bigfloat(x; precision=p))
    return rop
end

@inline function trunc!(rop::BigFloat, x::FloatLike)
    p = precision(rop)
    mpfr_trunc(rop, _as_bigfloat(x; precision=p))
    return rop
end

@inline function frac!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_frac(z, a, r), rop, x; rounding=rounding)
end

@inline function modf!(rop_int::BigFloat, rop_frac::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    p = precision(rop_int)
    mpfr_modf(rop_int, rop_frac, _as_bigfloat(x; precision=p), _rounding(rounding))
    return rop_int
end

@inline function fmodquo!(rop::BigFloat, quo::Ref{Clong}, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_fmodquo(rop, quo, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function remquo!(rop::BigFloat, quo::Ref{Clong}, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_remquo(rop, quo, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

# ============================================================================
# BigFloat: Advanced binary operations (fma, fms, hypot, reldiff, dim, min, max)
# ============================================================================

@inline function fma!(rop::BigFloat, x::FloatLike, y::FloatLike, z::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_fma(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), _as_bigfloat(z; precision=p), r)
    return rop
end

@inline function fms!(rop::BigFloat, x::FloatLike, y::FloatLike, z::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_fms(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), _as_bigfloat(z; precision=p), r)
    return rop
end

@inline function hypot!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_hypot(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function reldiff!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_reldiff(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function dim!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_dim(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function min!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_min(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

@inline function max!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_max(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

# ============================================================================
# BigFloat: Bit operations (mul_2exp!, div_2exp! - Unified, handles ui/si internally)
# ============================================================================

@inline function mul_2exp!(rop::BigFloat, x::FloatLike, n::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_mul_2ui(rop, _as_bigfloat(x; precision=precision(rop)), n, _rounding(rounding))
    return rop
end

@inline function mul_2exp!(rop::BigFloat, x::FloatLike, n::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_mul_2si(rop, _as_bigfloat(x; precision=precision(rop)), n, _rounding(rounding))
    return rop
end

@inline function div_2exp!(rop::BigFloat, x::FloatLike, n::T; rounding=MPFRRoundNearest) where {T <: NativeUnsignedInt}
    mpfr_div_2ui(rop, _as_bigfloat(x; precision=precision(rop)), n, _rounding(rounding))
    return rop
end

@inline function div_2exp!(rop::BigFloat, x::FloatLike, n::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_div_2si(rop, _as_bigfloat(x; precision=precision(rop)), n, _rounding(rounding))
    return rop
end

# ============================================================================
# BigFloat: Floating-point manipulation (nexttoward, nextabove, nextbelow, get_exp, set_exp, setsign, copysign)
# ============================================================================

@inline function nexttoward!(rop::BigFloat, y::FloatLike)
    mpfr_nexttoward(rop, _as_bigfloat(y; precision=precision(rop)))
    return rop
end

@inline function nextabove!(rop::BigFloat)
    mpfr_nextabove(rop)
    return rop
end

@inline function nextbelow!(rop::BigFloat)
    mpfr_nextbelow(rop)
    return rop
end

@inline function get_exp!(rop::BigFloat)
    return mpfr_get_exp(rop)
end

@inline function set_exp!(rop::BigFloat, e::T) where {T <: NativeSignedInt}
    return mpfr_set_exp(rop, e)
end

@inline function setsign!(rop::BigFloat, x::FloatLike, s::T; rounding=MPFRRoundNearest) where {T <: NativeSignedInt}
    mpfr_setsign(rop, _as_bigfloat(x; precision=precision(rop)), s, _rounding(rounding))
    return rop
end

@inline function copysign!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_copysign(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), r)
    return rop
end

# ============================================================================
# BigFloat: Mathematical constants (pi!, euler!, catalan!)
# ============================================================================

@inline function pi!(rop::BigFloat; rounding=MPFRRoundNearest)
    mpfr_const_pi(rop, _rounding(rounding))
    return rop
end

@inline function euler!(rop::BigFloat; rounding=MPFRRoundNearest)
    mpfr_const_euler(rop, _rounding(rounding))
    return rop
end

@inline function catalan!(rop::BigFloat; rounding=MPFRRoundNearest)
    mpfr_const_catalan(rop, _rounding(rounding))
    return rop
end
