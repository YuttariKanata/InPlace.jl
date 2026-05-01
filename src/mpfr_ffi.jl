# src/mpfr_ffi.jl
# Low-level MPFR bindings using MPFR_jll.

using MPFR_jll
using Base.MPFR:
    MPFRRoundingMode,
    MPFRRoundNearest,
    MPFRRoundToZero,
    MPFRRoundUp,
    MPFRRoundDown,
    MPFRRoundFromZero,
    MPFRRoundFaithful
using Base.GMP: CulongMax, ClongMax, CdoubleMax
using Libdl

const mpfr_prec_t = Clong
const mpfr_exp_t = Clong
const mpfr_flags_t = Cuint
const mpfr_rnd_t = MPFRRoundingMode

const MPFRFloat = BigFloat
const MPFRPtr = Ref{BigFloat}
const MPFRArrayPtr = Ptr{MPFRPtr}
const MPFRCharPtr = Ptr{UInt8}
const MPFREndPtr = Ref{MPFRCharPtr}
const MPFRCvoidPtr = Ptr{Cvoid}

const MPFRIntArg = Union{Int8, Int16, Int32, Int64}
const MPFRUIntArg = Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64}
const MPFRSizeArg = Csize_t == UInt32 ?
    Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32} :
    Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64}
const MPFRPrecArg = Union{Int8, Int16, Int32, Int64}
const MPFRExpArg = Union{Int8, Int16, Int32, Int64}
const MPFRFlagsArg = Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32}

# Julia argument type, C argument type.
const MPFR_TYPE_MAP = Dict(
    :F      => (BigFloat,                       :(Ref{BigFloat})),
    :Z      => (BigInt,                         :(Ref{BigInt})),
    :UI     => (MPFRUIntArg,                    :Culong),
    :SI     => (MPFRIntArg,                     :Clong),
    :UJ     => (MPFRUIntArg,                    :Cuintmax_t),
    :SJ     => (MPFRIntArg,                     :Cintmax_t),
    :I      => (MPFRIntArg,                     :Cint),
    :D      => (CdoubleMax,                     :Cdouble),
    :FLT    => (Union{Float16, Float32},        :Cfloat),
    :RND    => (MPFRRoundingMode,               :MPFRRoundingMode),
    :PREC   => (MPFRPrecArg,                    :mpfr_prec_t),
    :EXP    => (MPFRExpArg,                     :mpfr_exp_t),
    :SZ     => (MPFRSizeArg,                    :Csize_t),
    :FLAGS  => (MPFRFlagsArg,                   :mpfr_flags_t),
    :S      => (AbstractString,                 :Cstring),
    :PCHAR  => (Union{MPFRCharPtr, Ptr{Nothing}}, :(Ptr{UInt8})),
    :PPCHAR => (MPFREndPtr,                     :(Ptr{Ptr{UInt8}})),
    :PV     => (Ptr,                            :(Ptr{Cvoid})),
    :PINT   => (Ref{Cint},                      :(Ptr{Cint})),
    :PEXP   => (Ref{mpfr_exp_t},                :(Ptr{mpfr_exp_t})),
    :PCLONG => (Ref{Clong},                     :(Ptr{Clong})),
    :PFARR  => (MPFRArrayPtr,                   :(Ptr{Ref{BigFloat}})),
    :RAND   => (MPFRCvoidPtr,                   :(Ptr{Cvoid})),
)

const MPFR_RET_MAP = Dict(
    :Cvoid   => (:Nothing,       :Cvoid),
    :Cint    => (:Cint,          :Cint),
    :Cuint   => (:Cuint,         :Cuint),
    :Clong   => (:Clong,         :Clong),
    :Culong  => (:Culong,        :Culong),
    :Cintmax => (:Cintmax_t,     :Cintmax_t),
    :Cuintmax => (:Cuintmax_t,   :Cuintmax_t),
    :Cfloat  => (:Cfloat,        :Cfloat),
    :Cdouble => (:Cdouble,       :Cdouble),
    :Csize_t => (:Csize_t,       :Csize_t),
    :PREC    => (:mpfr_prec_t,   :mpfr_prec_t),
    :EXP     => (:mpfr_exp_t,    :mpfr_exp_t),
    :FLAGS   => (:mpfr_flags_t,  :mpfr_flags_t),
    :RND     => (:MPFRRoundingMode, :MPFRRoundingMode),
    :PCHAR   => (:MPFRCharPtr,   :(Ptr{UInt8})),
    :PV      => (:MPFRCvoidPtr,  :(Ptr{Cvoid})),
)

# These are deliberately not auto-bound. They either create/destroy MPFR
# objects owned by Julia, need unavailable C scalar types, need mpq_t/mpf_t,
# or are varargs interfaces.
const MPFR_UNSUPPORTED = Set([
    :mpfr_init,
    :mpfr_init2,
    :mpfr_inits,
    :mpfr_inits2,
    :mpfr_clear,
    :mpfr_clears,
    :mpfr_init_set,
    :mpfr_init_set_ui,
    :mpfr_init_set_si,
    :mpfr_init_set_d,
    :mpfr_init_set_ld,
    :mpfr_init_set_z,
    :mpfr_init_set_q,
    :mpfr_init_set_f,
    :mpfr_init_set_str,
    :mpfr_set_prec_raw,
    :mpfr_set_ld,
    :mpfr_get_ld,
    :mpfr_get_ld_2exp,
    :mpfr_set_float128,
    :mpfr_get_float128,
    :mpfr_set_decimal64,
    :mpfr_set_decimal128,
    :mpfr_get_decimal64,
    :mpfr_get_decimal128,
    :mpfr_set_q,
    :mpfr_set_f,
    :mpfr_get_q,
    :mpfr_get_f,
    :mpfr_add_q,
    :mpfr_sub_q,
    :mpfr_mul_q,
    :mpfr_div_q,
    :mpfr_cmp_q,
    :mpfr_cmp_f,
    :mpfr_fprintf,
    :mpfr_printf,
    :mpfr_sprintf,
    :mpfr_snprintf,
    :mpfr_asprintf,
    :mpfr_vfprintf,
    :mpfr_vprintf,
    :mpfr_vsprintf,
    :mpfr_vsnprintf,
    :mpfr_vasprintf,
])

# (return type, function name, argument key list)
const MPFR_DEFINITIONS = [
    # Initialization, precision and assignment for existing Julia BigFloat objects.
    (:Cvoid,   :mpfr_set_default_prec, [:PREC]),
    (:PREC,    :mpfr_get_default_prec, []),
    (:Cvoid,   :mpfr_set_default_rounding_mode, [:RND]),
    (:RND,     :mpfr_get_default_rounding_mode, []),
    (:Cvoid,   :mpfr_set_prec,         [:F, :PREC]),
    (:PREC,    :mpfr_get_prec,         [:F]),
    (:Cvoid,   :mpfr_prec_round,       [:F, :PREC, :RND]),
    (:Cint,    :mpfr_can_round,        [:F, :EXP, :RND, :RND, :PREC]),
    (:PREC,    :mpfr_min_prec,         [:F]),

    (:Cint,    :mpfr_set,              [:F, :F, :RND]),
    (:Cint,    :mpfr_set_ui,           [:F, :UI, :RND]),
    (:Cint,    :mpfr_set_si,           [:F, :SI, :RND]),
    (:Cint,    :mpfr_set_uj,           [:F, :UJ, :RND]),
    (:Cint,    :mpfr_set_sj,           [:F, :SJ, :RND]),
    (:Cint,    :mpfr_set_flt,          [:F, :FLT, :RND]),
    (:Cint,    :mpfr_set_d,            [:F, :D, :RND]),
    (:Cint,    :mpfr_set_z,            [:F, :Z, :RND]),
    (:Cint,    :mpfr_set_ui_2exp,      [:F, :UI, :EXP, :RND]),
    (:Cint,    :mpfr_set_si_2exp,      [:F, :SI, :EXP, :RND]),
    (:Cint,    :mpfr_set_uj_2exp,      [:F, :UJ, :SJ, :RND]),
    (:Cint,    :mpfr_set_sj_2exp,      [:F, :SJ, :SJ, :RND]),
    (:Cint,    :mpfr_set_z_2exp,       [:F, :Z, :EXP, :RND]),
    (:Cint,    :mpfr_set_str,          [:F, :S, :I, :RND]),
    (:Cint,    :mpfr_strtofr,          [:F, :S, :PPCHAR, :I, :RND]),
    (:Cvoid,   :mpfr_set_nan,          [:F]),
    (:Cvoid,   :mpfr_set_inf,          [:F, :I]),
    (:Cvoid,   :mpfr_set_zero,         [:F, :I]),
    (:Cvoid,   :mpfr_swap,             [:F, :F]),

    # Conversion functions.
    (:Cfloat,  :mpfr_get_flt,          [:F, :RND]),
    (:Cdouble, :mpfr_get_d,            [:F, :RND]),
    (:Clong,   :mpfr_get_si,           [:F, :RND]),
    (:Culong,  :mpfr_get_ui,           [:F, :RND]),
    (:Cintmax, :mpfr_get_sj,           [:F, :RND]),
    (:Cuintmax,:mpfr_get_uj,           [:F, :RND]),
    (:Cdouble, :mpfr_get_d_2exp,       [:PCLONG, :F, :RND]),
    (:Cint,    :mpfr_frexp,            [:PEXP, :F, :F, :RND]),
    (:EXP,     :mpfr_get_z_2exp,       [:Z, :F]),
    (:Cint,    :mpfr_get_z,            [:Z, :F, :RND]),
    (:Csize_t, :mpfr_get_str_ndigits,  [:I, :PREC]),
    (:PCHAR,   :mpfr_get_str,          [:PCHAR, :PEXP, :I, :SZ, :F, :RND]),
    (:Cvoid,   :mpfr_free_str,         [:PCHAR]),

    # Fits functions.
    (:Cint,    :mpfr_fits_ulong_p,     [:F, :RND]),
    (:Cint,    :mpfr_fits_slong_p,     [:F, :RND]),
    (:Cint,    :mpfr_fits_uint_p,      [:F, :RND]),
    (:Cint,    :mpfr_fits_sint_p,      [:F, :RND]),
    (:Cint,    :mpfr_fits_ushort_p,    [:F, :RND]),
    (:Cint,    :mpfr_fits_sshort_p,    [:F, :RND]),
    (:Cint,    :mpfr_fits_uintmax_p,   [:F, :RND]),
    (:Cint,    :mpfr_fits_intmax_p,    [:F, :RND]),

    # Arithmetic.
    (:Cint,    :mpfr_add,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_add_ui,           [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_add_si,           [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_add_d,            [:F, :F, :D, :RND]),
    (:Cint,    :mpfr_add_z,            [:F, :F, :Z, :RND]),
    (:Cint,    :mpfr_sub,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_ui_sub,           [:F, :UI, :F, :RND]),
    (:Cint,    :mpfr_sub_ui,           [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_si_sub,           [:F, :SI, :F, :RND]),
    (:Cint,    :mpfr_sub_si,           [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_d_sub,            [:F, :D, :F, :RND]),
    (:Cint,    :mpfr_sub_d,            [:F, :F, :D, :RND]),
    (:Cint,    :mpfr_z_sub,            [:F, :Z, :F, :RND]),
    (:Cint,    :mpfr_sub_z,            [:F, :F, :Z, :RND]),
    (:Cint,    :mpfr_mul,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_mul_ui,           [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_mul_si,           [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_mul_d,            [:F, :F, :D, :RND]),
    (:Cint,    :mpfr_mul_z,            [:F, :F, :Z, :RND]),
    (:Cint,    :mpfr_sqr,              [:F, :F, :RND]),
    (:Cint,    :mpfr_div,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_ui_div,           [:F, :UI, :F, :RND]),
    (:Cint,    :mpfr_div_ui,           [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_si_div,           [:F, :SI, :F, :RND]),
    (:Cint,    :mpfr_div_si,           [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_d_div,            [:F, :D, :F, :RND]),
    (:Cint,    :mpfr_div_d,            [:F, :F, :D, :RND]),
    (:Cint,    :mpfr_div_z,            [:F, :F, :Z, :RND]),
    (:Cint,    :mpfr_sqrt,             [:F, :F, :RND]),
    (:Cint,    :mpfr_sqrt_ui,          [:F, :UI, :RND]),
    (:Cint,    :mpfr_rec_sqrt,         [:F, :F, :RND]),
    (:Cint,    :mpfr_cbrt,             [:F, :F, :RND]),
    (:Cint,    :mpfr_rootn_ui,         [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_rootn_si,         [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_root,             [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_neg,              [:F, :F, :RND]),
    (:Cint,    :mpfr_abs,              [:F, :F, :RND]),
    (:Cint,    :mpfr_dim,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_mul_2ui,          [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_mul_2si,          [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_div_2ui,          [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_div_2si,          [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_mul_2exp,         [:F, :F, :EXP, :RND]),
    (:Cint,    :mpfr_div_2exp,         [:F, :F, :EXP, :RND]),
    (:Cint,    :mpfr_fac_ui,           [:F, :UI, :RND]),
    (:Cint,    :mpfr_fma,              [:F, :F, :F, :F, :RND]),
    (:Cint,    :mpfr_fms,              [:F, :F, :F, :F, :RND]),
    (:Cint,    :mpfr_fmma,             [:F, :F, :F, :F, :F, :RND]),
    (:Cint,    :mpfr_fmms,             [:F, :F, :F, :F, :F, :RND]),
    (:Cint,    :mpfr_hypot,            [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_reldiff,          [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_min,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_max,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_sum,              [:F, :PFARR, :UI, :RND]),
    (:Cint,    :mpfr_dot,              [:F, :PFARR, :PFARR, :UI, :RND]),

    # Comparison.
    (:Cint,    :mpfr_cmp,              [:F, :F]),
    (:Cint,    :mpfr_cmp_ui,           [:F, :UI]),
    (:Cint,    :mpfr_cmp_si,           [:F, :SI]),
    (:Cint,    :mpfr_cmp_d,            [:F, :D]),
    (:Cint,    :mpfr_cmp_z,            [:F, :Z]),
    (:Cint,    :mpfr_cmp_ui_2exp,      [:F, :UI, :EXP]),
    (:Cint,    :mpfr_cmp_si_2exp,      [:F, :SI, :EXP]),
    (:Cint,    :mpfr_cmpabs,           [:F, :F]),
    (:Cint,    :mpfr_cmpabs_ui,        [:F, :UI]),
    (:Cint,    :mpfr_nan_p,            [:F]),
    (:Cint,    :mpfr_inf_p,            [:F]),
    (:Cint,    :mpfr_number_p,         [:F]),
    (:Cint,    :mpfr_zero_p,           [:F]),
    (:Cint,    :mpfr_regular_p,        [:F]),
    (:Cint,    :mpfr_sgn,              [:F]),
    (:Cint,    :mpfr_greater_p,        [:F, :F]),
    (:Cint,    :mpfr_greaterequal_p,   [:F, :F]),
    (:Cint,    :mpfr_less_p,           [:F, :F]),
    (:Cint,    :mpfr_lessequal_p,      [:F, :F]),
    (:Cint,    :mpfr_equal_p,          [:F, :F]),
    (:Cint,    :mpfr_lessgreater_p,    [:F, :F]),
    (:Cint,    :mpfr_unordered_p,      [:F, :F]),
    (:Cint,    :mpfr_total_order_p,    [:F, :F]),

    # Elementary and special functions.
    (:Cint,    :mpfr_log,              [:F, :F, :RND]),
    (:Cint,    :mpfr_log_ui,           [:F, :UI, :RND]),
    (:Cint,    :mpfr_log2,             [:F, :F, :RND]),
    (:Cint,    :mpfr_log10,            [:F, :F, :RND]),
    (:Cint,    :mpfr_log1p,            [:F, :F, :RND]),
    (:Cint,    :mpfr_log2p1,           [:F, :F, :RND]),
    (:Cint,    :mpfr_log10p1,          [:F, :F, :RND]),
    (:Cint,    :mpfr_exp,              [:F, :F, :RND]),
    (:Cint,    :mpfr_exp2,             [:F, :F, :RND]),
    (:Cint,    :mpfr_exp10,            [:F, :F, :RND]),
    (:Cint,    :mpfr_expm1,            [:F, :F, :RND]),
    (:Cint,    :mpfr_exp2m1,           [:F, :F, :RND]),
    (:Cint,    :mpfr_exp10m1,          [:F, :F, :RND]),
    (:Cint,    :mpfr_pow,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_powr,             [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_pow_ui,           [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_pow_si,           [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_pow_uj,           [:F, :F, :UJ, :RND]),
    (:Cint,    :mpfr_pow_sj,           [:F, :F, :SJ, :RND]),
    (:Cint,    :mpfr_pown,             [:F, :F, :SJ, :RND]),
    (:Cint,    :mpfr_pow_z,            [:F, :F, :Z, :RND]),
    (:Cint,    :mpfr_ui_pow_ui,        [:F, :UI, :UI, :RND]),
    (:Cint,    :mpfr_ui_pow,           [:F, :UI, :F, :RND]),
    (:Cint,    :mpfr_compound_si,      [:F, :F, :SI, :RND]),
    (:Cint,    :mpfr_cos,              [:F, :F, :RND]),
    (:Cint,    :mpfr_sin,              [:F, :F, :RND]),
    (:Cint,    :mpfr_tan,              [:F, :F, :RND]),
    (:Cint,    :mpfr_cosu,             [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_sinu,             [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_tanu,             [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_cospi,            [:F, :F, :RND]),
    (:Cint,    :mpfr_sinpi,            [:F, :F, :RND]),
    (:Cint,    :mpfr_tanpi,            [:F, :F, :RND]),
    (:Cint,    :mpfr_sin_cos,          [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_sec,              [:F, :F, :RND]),
    (:Cint,    :mpfr_csc,              [:F, :F, :RND]),
    (:Cint,    :mpfr_cot,              [:F, :F, :RND]),
    (:Cint,    :mpfr_acos,             [:F, :F, :RND]),
    (:Cint,    :mpfr_asin,             [:F, :F, :RND]),
    (:Cint,    :mpfr_atan,             [:F, :F, :RND]),
    (:Cint,    :mpfr_acosu,            [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_asinu,            [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_atanu,            [:F, :F, :UI, :RND]),
    (:Cint,    :mpfr_acospi,           [:F, :F, :RND]),
    (:Cint,    :mpfr_asinpi,           [:F, :F, :RND]),
    (:Cint,    :mpfr_atanpi,           [:F, :F, :RND]),
    (:Cint,    :mpfr_atan2,            [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_atan2u,           [:F, :F, :F, :UI, :RND]),
    (:Cint,    :mpfr_atan2pi,          [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_cosh,             [:F, :F, :RND]),
    (:Cint,    :mpfr_sinh,             [:F, :F, :RND]),
    (:Cint,    :mpfr_tanh,             [:F, :F, :RND]),
    (:Cint,    :mpfr_sinh_cosh,        [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_sech,             [:F, :F, :RND]),
    (:Cint,    :mpfr_csch,             [:F, :F, :RND]),
    (:Cint,    :mpfr_coth,             [:F, :F, :RND]),
    (:Cint,    :mpfr_acosh,            [:F, :F, :RND]),
    (:Cint,    :mpfr_asinh,            [:F, :F, :RND]),
    (:Cint,    :mpfr_atanh,            [:F, :F, :RND]),
    (:Cint,    :mpfr_eint,             [:F, :F, :RND]),
    (:Cint,    :mpfr_li2,              [:F, :F, :RND]),
    (:Cint,    :mpfr_gamma,            [:F, :F, :RND]),
    (:Cint,    :mpfr_gamma_inc,        [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_lngamma,          [:F, :F, :RND]),
    (:Cint,    :mpfr_lgamma,           [:F, :PINT, :F, :RND]),
    (:Cint,    :mpfr_digamma,          [:F, :F, :RND]),
    (:Cint,    :mpfr_beta,             [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_zeta,             [:F, :F, :RND]),
    (:Cint,    :mpfr_zeta_ui,          [:F, :UI, :RND]),
    (:Cint,    :mpfr_erf,              [:F, :F, :RND]),
    (:Cint,    :mpfr_erfc,             [:F, :F, :RND]),
    (:Cint,    :mpfr_j0,               [:F, :F, :RND]),
    (:Cint,    :mpfr_j1,               [:F, :F, :RND]),
    (:Cint,    :mpfr_jn,               [:F, :SI, :F, :RND]),
    (:Cint,    :mpfr_y0,               [:F, :F, :RND]),
    (:Cint,    :mpfr_y1,               [:F, :F, :RND]),
    (:Cint,    :mpfr_yn,               [:F, :SI, :F, :RND]),
    (:Cint,    :mpfr_agm,              [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_ai,               [:F, :F, :RND]),
    (:Cint,    :mpfr_const_log2,       [:F, :RND]),
    (:Cint,    :mpfr_const_pi,         [:F, :RND]),
    (:Cint,    :mpfr_const_euler,      [:F, :RND]),
    (:Cint,    :mpfr_const_catalan,    [:F, :RND]),

    # Integer and remainder related functions.
    (:Cint,    :mpfr_rint,             [:F, :F, :RND]),
    (:Cint,    :mpfr_ceil,             [:F, :F]),
    (:Cint,    :mpfr_floor,            [:F, :F]),
    (:Cint,    :mpfr_round,            [:F, :F]),
    (:Cint,    :mpfr_roundeven,        [:F, :F]),
    (:Cint,    :mpfr_trunc,            [:F, :F]),
    (:Cint,    :mpfr_rint_ceil,        [:F, :F, :RND]),
    (:Cint,    :mpfr_rint_floor,       [:F, :F, :RND]),
    (:Cint,    :mpfr_rint_round,       [:F, :F, :RND]),
    (:Cint,    :mpfr_rint_roundeven,   [:F, :F, :RND]),
    (:Cint,    :mpfr_rint_trunc,       [:F, :F, :RND]),
    (:Cint,    :mpfr_frac,             [:F, :F, :RND]),
    (:Cint,    :mpfr_modf,             [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_fmod,             [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_fmodquo,          [:F, :PCLONG, :F, :F, :RND]),
    (:Cint,    :mpfr_remainder,        [:F, :F, :F, :RND]),
    (:Cint,    :mpfr_remquo,           [:F, :PCLONG, :F, :F, :RND]),
    (:Cint,    :mpfr_integer_p,        [:F]),

    # Miscellaneous floating-point manipulation.
    (:Cint,    :mpfr_nexttoward,       [:F, :F]),
    (:Cvoid,   :mpfr_nextabove,        [:F]),
    (:Cvoid,   :mpfr_nextbelow,        [:F]),
    (:EXP,     :mpfr_get_exp,          [:F]),
    (:Cint,    :mpfr_set_exp,          [:F, :EXP]),
    (:Cint,    :mpfr_signbit,          [:F]),
    (:Cint,    :mpfr_setsign,          [:F, :F, :I, :RND]),
    (:Cint,    :mpfr_copysign,         [:F, :F, :F, :RND]),
    # Input and output through FILE* represented as Ptr{Cvoid}.
    (:Csize_t, :mpfr_out_str,          [:PV, :I, :SZ, :F, :RND]),
    (:Csize_t, :mpfr_inp_str,          [:F, :PV, :I, :RND]),
    (:Cint,    :mpfr_fpif_export,      [:PV, :F]),
    (:Cint,    :mpfr_fpif_import,      [:F, :PV]),
    (:Cvoid,   :mpfr_dump,             [:F]),

    # Random generation.
    (:Cint,    :mpfr_urandomb,         [:F, :RAND]),
    (:Cint,    :mpfr_urandom,          [:F, :RAND, :RND]),
    (:Cint,    :mpfr_nrandom,          [:F, :RAND, :RND]),
    (:Cint,    :mpfr_erandom,          [:F, :RAND, :RND]),
    (:Cint,    :mpfr_grandom,          [:F, :F, :RAND, :RND]),

    # Exception flags and exponent range.
    (:Cvoid,   :mpfr_clear_flags,      []),
    (:Cvoid,   :mpfr_clear_underflow,  []),
    (:Cvoid,   :mpfr_clear_overflow,   []),
    (:Cvoid,   :mpfr_clear_divby0,     []),
    (:Cvoid,   :mpfr_clear_nanflag,    []),
    (:Cvoid,   :mpfr_clear_inexflag,   []),
    (:Cvoid,   :mpfr_clear_erangeflag, []),
    (:Cint,    :mpfr_underflow_p,      []),
    (:Cint,    :mpfr_overflow_p,       []),
    (:Cint,    :mpfr_divby0_p,         []),
    (:Cint,    :mpfr_nanflag_p,        []),
    (:Cint,    :mpfr_inexflag_p,       []),
    (:Cint,    :mpfr_erangeflag_p,     []),
    (:FLAGS,   :mpfr_flags_save,       []),
    (:Cvoid,   :mpfr_flags_restore,    [:FLAGS, :I]),
    (:Cvoid,   :mpfr_flags_clear,      [:FLAGS]),
    (:Cvoid,   :mpfr_flags_set,        [:FLAGS]),
    (:Cint,    :mpfr_flags_test,       [:FLAGS]),
    (:Cvoid,   :mpfr_set_underflow,    []),
    (:Cvoid,   :mpfr_set_overflow,     []),
    (:Cvoid,   :mpfr_set_divby0,       []),
    (:Cvoid,   :mpfr_set_nanflag,      []),
    (:Cvoid,   :mpfr_set_inexflag,     []),
    (:Cvoid,   :mpfr_set_erangeflag,   []),
    (:EXP,     :mpfr_get_emin,         []),
    (:EXP,     :mpfr_get_emax,         []),
    (:Cint,    :mpfr_set_emin,         [:EXP]),
    (:Cint,    :mpfr_set_emax,         [:EXP]),
    (:EXP,     :mpfr_get_emin_min,     []),
    (:EXP,     :mpfr_get_emin_max,     []),
    (:EXP,     :mpfr_get_emax_min,     []),
    (:EXP,     :mpfr_get_emax_max,     []),
    (:Cint,    :mpfr_check_range,      [:F, :I, :RND]),
    (:Cint,    :mpfr_subnormalize,     [:F, :I, :RND]),

    # Version information and build metadata.
    (:PCHAR,   :mpfr_get_version,      []),
    (:PCHAR,   :mpfr_get_patches,      []),
    (:Cint,    :mpfr_buildopt_tls_p,   []),
    (:Cint,    :mpfr_buildopt_float128_p, []),
    (:Cint,    :mpfr_buildopt_decimal_p, []),
    (:Cint,    :mpfr_buildopt_gmpinternals_p, []),
    (:PCHAR,   :mpfr_buildopt_sharedcache_p, []),
    (:PCHAR,   :mpfr_buildopt_tune_case, []),
    (:Cvoid,   :mpfr_free_cache,       []),
    (:Cvoid,   :mpfr_free_cache2,      [:I]),
    (:Cvoid,   :mpfr_free_pool,        []),
    (:Cvoid,   :mpfr_mp_memory_cleanup, []),

    # Custom interface. These are low-level and exposed as raw pointers/integers.
    (:Csize_t, :mpfr_custom_get_size,  [:PREC]),
    (:Cvoid,   :mpfr_custom_init,      [:PV, :PREC]),
    (:Cvoid,   :mpfr_custom_init_set,  [:F, :I, :EXP, :PREC, :PV]),
    (:PV,      :mpfr_custom_get_significand, [:F]),
    (:EXP,     :mpfr_custom_get_exp,   [:F]),
    (:Cint,    :mpfr_custom_get_kind,  [:F]),
    (:Cvoid,   :mpfr_custom_move,      [:F, :PV]),
]

const MPFR_SYMBOL_OVERRIDES = Dict(
    :mpfr_get_sj => :__gmpfr_mpfr_get_sj,
    :mpfr_get_uj => :__gmpfr_mpfr_get_uj,
)

function mpfr_symbol(fname::Symbol)
    get(MPFR_SYMBOL_OVERRIDES, fname, fname)
end

function has_mpfr_symbol(fname::Symbol)
    ptr = dlsym_e(dlopen(libmpfr), mpfr_symbol(fname))
    return ptr != C_NULL
end

mpfr_arg(::Val, x) = x
mpfr_arg(::Val{:UI}, x) = Culong(x)
mpfr_arg(::Val{:SI}, x) = Clong(x)
mpfr_arg(::Val{:UJ}, x) = Cuintmax_t(x)
mpfr_arg(::Val{:SJ}, x) = Cintmax_t(x)
mpfr_arg(::Val{:I}, x) = Cint(x)
mpfr_arg(::Val{:D}, x) = Cdouble(x)
mpfr_arg(::Val{:FLT}, x) = Cfloat(x)
mpfr_arg(::Val{:PREC}, x) = mpfr_prec_t(x)
mpfr_arg(::Val{:EXP}, x) = mpfr_exp_t(x)
mpfr_arg(::Val{:SZ}, x) = Csize_t(x)
mpfr_arg(::Val{:FLAGS}, x) = mpfr_flags_t(x)

const MPFR_MISSING = Symbol[]

function mpfr_ffi_defining_functions()
    empty!(MPFR_MISSING)
    for (ret_type, fname, arg_keys) in MPFR_DEFINITIONS
        if !has_mpfr_symbol(fname)
            push!(MPFR_MISSING, fname)
            continue
        end

        j_types = [MPFR_TYPE_MAP[k][1] for k in arg_keys]
        c_types = [MPFR_TYPE_MAP[k][2] for k in arg_keys]
        args = [Symbol(:a, i) for i in eachindex(arg_keys)]
        c_args = [:($(mpfr_arg)(Val{$(QuoteNode(arg_keys[i]))}(), $(args[i]))) for i in eachindex(args)]
        sig = [Expr(:(::), args[i], j_types[i]) for i in eachindex(args)]
        j_ret = MPFR_RET_MAP[ret_type][1]
        c_ret = MPFR_RET_MAP[ret_type][2]
        target = (mpfr_symbol(fname), libmpfr)

        @eval begin
            export $fname
            function $fname($(sig...))::$j_ret
                return ccall($target, $c_ret, ($(c_types...),), $(c_args...))
            end
        end
    end
    return MPFR_MISSING
end

mpfr_ffi_defining_functions()

export MPFRRoundingMode,
    MPFRRoundNearest,
    MPFRRoundToZero,
    MPFRRoundUp,
    MPFRRoundDown,
    MPFRRoundFromZero,
    MPFRRoundFaithful,
    mpfr_prec_t,
    mpfr_exp_t,
    mpfr_flags_t,
    MPFR_MISSING
