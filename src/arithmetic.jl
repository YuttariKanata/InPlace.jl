# User-facing mutating arithmetic helpers.

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
    sin!,
    cos!,
    tan!,
    exp!,
    log!,
    min!,
    max!,
    pi!,
    euler!,
    catalan!

const IntegerLike = Union{Integer, BigInt}
const FloatLike = Union{AbstractFloat, Integer, BigInt, BigFloat}

_rounding(r::MPFRRoundingMode) = r
_rounding(r::RoundingMode) = convert(MPFRRoundingMode, r)

_as_bigint(x::BigInt) = x
_as_bigint(x::Integer) = BigInt(x)

function _checked_nonnegative(::Val{name}, x::Integer) where {name}
    x < 0 && throw(ArgumentError("$(name) requires a non-negative integer argument"))
    return x
end

function _as_bigfloat(x::BigFloat)
    return x
end

function _as_bigfloat(x::FloatLike; precision::Integer=precision(BigFloat))
    return BigFloat(x; precision=precision)
end

function _with_bigfloat_arg(f, rop::BigFloat, x; rounding=MPFRRoundNearest)
    return f(rop, _as_bigfloat(x; precision=precision(rop)), _rounding(rounding))
end

# BigInt assignment and arithmetic.

function set!(rop::BigInt, x::BigInt)
    mpz_set(rop, x)
    return rop
end

function set!(rop::BigInt, x::Integer)
    if x < 0
        mpz_set_si(rop, x)
    else
        mpz_set_ui(rop, x)
    end
    return rop
end

function set!(rop::BigInt, s::AbstractString; base::Integer=0)
    rc = mpz_set_str(rop, s, base)
    rc == 0 || throw(ArgumentError("invalid integer string for base $base: $s"))
    return rop
end

function swap!(x::BigInt, y::BigInt)
    mpz_swap(x, y)
    return x
end

function add!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_add(rop, x, y)
    return rop
end

function add!(rop::BigInt, x::BigInt, y::Integer)
    if y < 0
        mpz_sub_ui(rop, x, -y)
    else
        mpz_add_ui(rop, x, y)
    end
    return rop
end

add!(rop::BigInt, x::Integer, y::BigInt) = add!(rop, y, x)

function add!(rop::BigInt, x::Integer, y::Integer)
    return set!(rop, _as_bigint(x) + _as_bigint(y))
end

function sub!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_sub(rop, x, y)
    return rop
end

function sub!(rop::BigInt, x::BigInt, y::Integer)
    if y < 0
        mpz_add_ui(rop, x, -y)
    else
        mpz_sub_ui(rop, x, y)
    end
    return rop
end

function sub!(rop::BigInt, x::Integer, y::BigInt)
    if x < 0
        mpz_add_ui(rop, y, -x)
        mpz_neg(rop, rop)
    else
        mpz_ui_sub(rop, x, y)
    end
    return rop
end

function sub!(rop::BigInt, x::Integer, y::Integer)
    return set!(rop, _as_bigint(x) - _as_bigint(y))
end

function mul!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_mul(rop, x, y)
    return rop
end

function mul!(rop::BigInt, x::BigInt, y::Integer)
    if y < 0
        mpz_mul_ui(rop, x, -y)
        mpz_neg(rop, rop)
    else
        mpz_mul_ui(rop, x, y)
    end
    return rop
end

mul!(rop::BigInt, x::Integer, y::BigInt) = mul!(rop, y, x)

function mul!(rop::BigInt, x::Integer, y::Integer)
    return set!(rop, _as_bigint(x) * _as_bigint(y))
end

function div!(rop::BigInt, x::BigInt, y::BigInt; rounding::Symbol=:trunc)
    rounding === :trunc && mpz_tdiv_q(rop, x, y)
    rounding === :floor && mpz_fdiv_q(rop, x, y)
    rounding === :ceil && mpz_cdiv_q(rop, x, y)
    rounding in (:trunc, :floor, :ceil) || throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    return rop
end

function div!(rop::BigInt, x::IntegerLike, y::IntegerLike; rounding::Symbol=:trunc)
    return div!(rop, _as_bigint(x), _as_bigint(y); rounding=rounding)
end

function rem!(rop::BigInt, x::BigInt, y::BigInt; rounding::Symbol=:trunc)
    rounding === :trunc && mpz_tdiv_r(rop, x, y)
    rounding === :floor && mpz_fdiv_r(rop, x, y)
    rounding === :ceil && mpz_cdiv_r(rop, x, y)
    rounding in (:trunc, :floor, :ceil) || throw(ArgumentError("rounding must be :trunc, :floor, or :ceil"))
    return rop
end

function rem!(rop::BigInt, x::IntegerLike, y::IntegerLike; rounding::Symbol=:trunc)
    return rem!(rop, _as_bigint(x), _as_bigint(y); rounding=rounding)
end

function mod!(rop::BigInt, x::BigInt, y::BigInt)
    mpz_mod(rop, x, y)
    return rop
end

mod!(rop::BigInt, x::IntegerLike, y::IntegerLike) = mod!(rop, _as_bigint(x), _as_bigint(y))

function pow!(rop::BigInt, x::BigInt, y::Integer)
    _checked_nonnegative(Val(:pow!), y)
    mpz_pow_ui(rop, x, y)
    return rop
end

function pow!(rop::BigInt, x::Integer, y::Integer)
    _checked_nonnegative(Val(:pow!), y)
    if x < 0
        mpz_ui_pow_ui(rop, -x, y)
        isodd(y) && mpz_neg(rop, rop)
    else
        mpz_ui_pow_ui(rop, x, y)
    end
    return rop
end

function sqrt!(rop::BigInt, x::IntegerLike)
    mpz_sqrt(rop, _as_bigint(x))
    return rop
end

function neg!(rop::BigInt, x::IntegerLike)
    mpz_neg(rop, _as_bigint(x))
    return rop
end

function abs!(rop::BigInt, x::IntegerLike)
    mpz_abs(rop, _as_bigint(x))
    return rop
end

function gcd!(rop::BigInt, x::IntegerLike, y::IntegerLike)
    mpz_gcd(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

function lcm!(rop::BigInt, x::IntegerLike, y::IntegerLike)
    mpz_lcm(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

function and!(rop::BigInt, x::IntegerLike, y::IntegerLike)
    mpz_and(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

function or!(rop::BigInt, x::IntegerLike, y::IntegerLike)
    mpz_ior(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

function xor!(rop::BigInt, x::IntegerLike, y::IntegerLike)
    mpz_xor(rop, _as_bigint(x), _as_bigint(y))
    return rop
end

function com!(rop::BigInt, x::IntegerLike)
    mpz_com(rop, _as_bigint(x))
    return rop
end

# BigFloat assignment and arithmetic.

function set!(rop::BigFloat, x::BigFloat; rounding=MPFRRoundNearest)
    mpfr_set(rop, x, _rounding(rounding))
    return rop
end

function set!(rop::BigFloat, x::Integer; rounding=MPFRRoundNearest)
    if x isa BigInt
        mpfr_set_z(rop, x, _rounding(rounding))
    elseif x < 0
        mpfr_set_si(rop, x, _rounding(rounding))
    else
        mpfr_set_ui(rop, x, _rounding(rounding))
    end
    return rop
end

function set!(rop::BigFloat, x::AbstractFloat; rounding=MPFRRoundNearest)
    mpfr_set_d(rop, x, _rounding(rounding))
    return rop
end

function set!(rop::BigFloat, s::AbstractString; base::Integer=0, rounding=MPFRRoundNearest)
    rc = mpfr_set_str(rop, s, base, _rounding(rounding))
    rc == 0 || throw(ArgumentError("invalid floating-point string for base $base: $s"))
    return rop
end

function swap!(x::BigFloat, y::BigFloat)
    mpfr_swap(x, y)
    return x
end

function add!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_add(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function sub!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_sub(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function mul!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_mul(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function div!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_div(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function square!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sqr(z, a, r), rop, x; rounding=rounding)
end

function sqrt!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_sqrt(z, a, r), rop, x; rounding=rounding)
end

function neg!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_neg(z, a, r), rop, x; rounding=rounding)
end

function abs!(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
    return _with_bigfloat_arg((z, a, r) -> mpfr_abs(z, a, r), rop, x; rounding=rounding)
end

function pow!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_pow(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function fma!(rop::BigFloat, x::FloatLike, y::FloatLike, z::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    p = precision(rop)
    mpfr_fma(rop, _as_bigfloat(x; precision=p), _as_bigfloat(y; precision=p), _as_bigfloat(z; precision=p), r)
    return rop
end

function rem!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_remainder(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function mod!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_fmod(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function min!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_min(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

function max!(rop::BigFloat, x::FloatLike, y::FloatLike; rounding=MPFRRoundNearest)
    r = _rounding(rounding)
    mpfr_max(rop, _as_bigfloat(x; precision=precision(rop)), _as_bigfloat(y; precision=precision(rop)), r)
    return rop
end

for (jname, cname) in (
    (:sin!, :mpfr_sin),
    (:cos!, :mpfr_cos),
    (:tan!, :mpfr_tan),
    (:exp!, :mpfr_exp),
    (:log!, :mpfr_log),
)
    @eval begin
        function $jname(rop::BigFloat, x::FloatLike; rounding=MPFRRoundNearest)
            return _with_bigfloat_arg((z, a, r) -> $cname(z, a, r), rop, x; rounding=rounding)
        end
    end
end

function pi!(rop::BigFloat; rounding=MPFRRoundNearest)
    mpfr_const_pi(rop, _rounding(rounding))
    return rop
end

function euler!(rop::BigFloat; rounding=MPFRRoundNearest)
    mpfr_const_euler(rop, _rounding(rounding))
    return rop
end

function catalan!(rop::BigFloat; rounding=MPFRRoundNearest)
    mpfr_const_catalan(rop, _rounding(rounding))
    return rop
end
