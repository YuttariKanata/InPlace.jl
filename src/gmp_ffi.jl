# src/gmp_ffi.jl
# Low-level GMP mpz bindings using GMP_jll.

using GMP_jll
using Base.GMP.MPZ: mpz_t, bitcnt_t
using Base.GMP: CulongMax, ClongMax, CdoubleMax, Limb
using Libdl

const Cbitcnt = bitcnt_t
const mp_limb_t = Limb
const mp_size_t = BigInt.types[2]

const CintMax = Cint == Int32 ? Union{Int8, Int16, Int32} : Union{Int8, Int16, Int32, Int64}
const CulongArg = CulongMax
const ClongArg = ClongMax
const CbitcntMax = bitcnt_t == UInt32 ? Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32} : Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64}
const Csize_tMax = Csize_t == UInt32 ? Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32} : Union{Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64}
const Cmp_size_tMax = CintMax
const Cmp_limb_tMax = mp_limb_t == UInt32 ? Union{UInt8, UInt16, UInt32} : Union{UInt8, UInt16, UInt32, UInt64}

const CcharPtr = Ptr{UInt8}
const CvoidPtr = Ptr{Cvoid}
const LimbPtr = Ptr{mp_limb_t}

# Julia argument type, C argument type.
const TYPE_MAP = Dict(
    :M     => (BigInt,                         :mpz_t),
    :UI    => (CulongArg,                      :Culong),
    :SI    => (ClongArg,                       :Clong),
    :I     => (CintMax,                        :Cint),
    :D     => (CdoubleMax,                     :Cdouble),
    :BC    => (CbitcntMax,                     :bitcnt_t),
    :S     => (AbstractString,                 :Cstring),
    :PCHAR => (Union{CcharPtr, Ptr{Nothing}},  :(Ptr{UInt8})),
    :PV    => (Ptr,                            :(Ptr{Cvoid})),
    :PSI   => (Ref{Clong},                     :(Ptr{Clong})),
    :PUI   => (Ref{Culong},                    :(Ptr{Culong})),
    :PSZ   => (Ref{Csize_t},                   :(Ptr{Csize_t})),
    :SZ    => (Csize_tMax,                     :Csize_t),
    :ML    => (Cmp_limb_tMax,                  :mp_limb_t),
    :PML   => (LimbPtr,                        :(Ptr{mp_limb_t})),
    :MS    => (Cmp_size_tMax,                  :mp_size_t),
    :RAND  => (CvoidPtr,                       :(Ptr{Cvoid})),
)

const RET_MAP = Dict(
    :Cint      => (:Cint,      :Cint),
    :Cvoid     => (:Nothing,   :Cvoid),
    :Culong    => (:Culong,    :Culong),
    :Clong     => (:Clong,     :Clong),
    :Cdouble   => (:Cdouble,   :Cdouble),
    :Cbitcnt   => (:Cbitcnt,   :bitcnt_t),
    :Csize_t   => (:Csize_t,   :Csize_t),
    :mp_limb_t => (:mp_limb_t, :mp_limb_t),
    :PCHAR     => (:CcharPtr,  :(Ptr{UInt8})),
    :PV        => (:CvoidPtr,  :(Ptr{Cvoid})),
    :PML       => (:LimbPtr,   :(Ptr{mp_limb_t})),
)

# Functions that allocate, clear, or initialize an mpz_t are intentionally not
# bound for BigInt. Julia owns BigInt initialization and finalization.
const GMP_MPZ_UNSUPPORTED = Set([
    :mpz_init,
    :mpz_init2,
    :mpz_init_set,
    :mpz_init_set_d,
    :mpz_init_set_si,
    :mpz_init_set_str,
    :mpz_init_set_ui,
    :mpz_inits,
    :mpz_clear,
    :mpz_clears,
    :mpz_array_init,
    :mpz_roinit_n,
    # These require GMP mpq_t/mpf_t bindings, not Julia BigInt.
    :mpz_set_q,
    :mpz_set_f,
])

# (return type, function name, argument key list)
const GMP_DEFINITIONS = [
    # 5.1 Initialization Functions, safe for existing Julia BigInt objects.
    (:Cvoid,     :mpz_realloc2,             [:M,     :BC                    ]),

    # 5.2 Assignment Functions
    (:Cvoid,     :mpz_set,                  [:M,     :M                     ]),
    (:Cvoid,     :mpz_set_ui,               [:M,     :UI                    ]),
    (:Cvoid,     :mpz_set_si,               [:M,     :SI                    ]),
    (:Cvoid,     :mpz_set_d,                [:M,     :D                     ]),
    (:Cint,      :mpz_set_str,              [:M,     :S,     :I             ]),
    (:Cvoid,     :mpz_swap,                 [:M,     :M                     ]),

    # 5.4 Conversion Functions
    (:Culong,    :mpz_get_ui,               [:M                            ]),
    (:Clong,     :mpz_get_si,               [:M                            ]),
    (:Cdouble,   :mpz_get_d,                [:M                            ]),
    (:Cdouble,   :mpz_get_d_2exp,           [:PSI,   :M                     ]),
    (:PCHAR,     :mpz_get_str,              [:PCHAR, :I,     :M             ]),

    # 5.5 Arithmetic Functions
    (:Cvoid,     :mpz_add,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_add_ui,               [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_sub,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_sub_ui,               [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_ui_sub,               [:M,     :UI,    :M             ]),
    (:Cvoid,     :mpz_mul,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_mul_si,               [:M,     :M,     :SI            ]),
    (:Cvoid,     :mpz_mul_ui,               [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_addmul,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_addmul_ui,            [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_submul,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_submul_ui,            [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_mul_2exp,             [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_neg,                  [:M,     :M                     ]),
    (:Cvoid,     :mpz_abs,                  [:M,     :M                     ]),

    # 5.6 Division Functions
    (:Cvoid,     :mpz_cdiv_q,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_cdiv_r,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_cdiv_qr,              [:M,     :M,     :M,     :M     ]),
    (:Culong,    :mpz_cdiv_q_ui,            [:M,     :M,     :UI            ]),
    (:Culong,    :mpz_cdiv_r_ui,            [:M,     :M,     :UI            ]),
    (:Culong,    :mpz_cdiv_qr_ui,           [:M,     :M,     :M,     :UI    ]),
    (:Culong,    :mpz_cdiv_ui,              [:M,     :UI                    ]),
    (:Cvoid,     :mpz_cdiv_q_2exp,          [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_cdiv_r_2exp,          [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_fdiv_q,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_fdiv_r,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_fdiv_qr,              [:M,     :M,     :M,     :M     ]),
    (:Culong,    :mpz_fdiv_q_ui,            [:M,     :M,     :UI            ]),
    (:Culong,    :mpz_fdiv_r_ui,            [:M,     :M,     :UI            ]),
    (:Culong,    :mpz_fdiv_qr_ui,           [:M,     :M,     :M,     :UI    ]),
    (:Culong,    :mpz_fdiv_ui,              [:M,     :UI                    ]),
    (:Cvoid,     :mpz_fdiv_q_2exp,          [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_fdiv_r_2exp,          [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_tdiv_q,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_tdiv_r,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_tdiv_qr,              [:M,     :M,     :M,     :M     ]),
    (:Culong,    :mpz_tdiv_q_ui,            [:M,     :M,     :UI            ]),
    (:Culong,    :mpz_tdiv_r_ui,            [:M,     :M,     :UI            ]),
    (:Culong,    :mpz_tdiv_qr_ui,           [:M,     :M,     :M,     :UI    ]),
    (:Culong,    :mpz_tdiv_ui,              [:M,     :UI                    ]),
    (:Cvoid,     :mpz_tdiv_q_2exp,          [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_tdiv_r_2exp,          [:M,     :M,     :BC            ]),
    (:Cvoid,     :mpz_mod,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_divexact,             [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_divexact_ui,          [:M,     :M,     :UI            ]),
    (:Cint,      :mpz_divisible_p,          [:M,     :M                     ]),
    (:Cint,      :mpz_divisible_ui_p,       [:M,     :UI                    ]),
    (:Cint,      :mpz_divisible_2exp_p,     [:M,     :BC                    ]),
    (:Cint,      :mpz_congruent_p,          [:M,     :M,     :M             ]),
    (:Cint,      :mpz_congruent_ui_p,       [:M,     :UI,    :UI            ]),
    (:Cint,      :mpz_congruent_2exp_p,     [:M,     :M,     :BC            ]),

    # 5.7 Exponentiation Functions
    (:Cvoid,     :mpz_powm,                 [:M,     :M,     :M,     :M     ]),
    (:Cvoid,     :mpz_powm_ui,              [:M,     :M,     :UI,    :M     ]),
    (:Cvoid,     :mpz_powm_sec,             [:M,     :M,     :M,     :M     ]),
    (:Cvoid,     :mpz_pow_ui,               [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_ui_pow_ui,            [:M,     :UI,    :UI            ]),

    # 5.8 Root Extraction Functions
    (:Cint,      :mpz_root,                 [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_rootrem,              [:M,     :M,     :M,     :UI    ]),
    (:Cvoid,     :mpz_sqrt,                 [:M,     :M                     ]),
    (:Cvoid,     :mpz_sqrtrem,              [:M,     :M,     :M             ]),
    (:Cint,      :mpz_perfect_power_p,      [:M                            ]),
    (:Cint,      :mpz_perfect_square_p,     [:M                            ]),

    # 5.9 Number Theoretic Functions
    (:Cint,      :mpz_probab_prime_p,       [:M,     :I                     ]),
    (:Cvoid,     :mpz_nextprime,            [:M,     :M                     ]),
    (:Cint,      :mpz_prevprime,            [:M,     :M                     ]),
    (:Cvoid,     :mpz_gcd,                  [:M,     :M,     :M             ]),
    (:Culong,    :mpz_gcd_ui,               [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_gcdext,               [:M,     :M,     :M,     :M,     :M]),
    (:Cvoid,     :mpz_lcm,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_lcm_ui,               [:M,     :M,     :UI            ]),
    (:Cint,      :mpz_invert,               [:M,     :M,     :M             ]),
    (:Cint,      :mpz_jacobi,               [:M,     :M                     ]),
    (:Cint,      :mpz_legendre,             [:M,     :M                     ]),
    (:Cint,      :mpz_kronecker_si,         [:M,     :SI                    ]),
    (:Cint,      :mpz_kronecker_ui,         [:M,     :UI                    ]),
    (:Cint,      :mpz_si_kronecker,         [:SI,    :M                     ]),
    (:Cint,      :mpz_ui_kronecker,         [:UI,    :M                     ]),
    (:Cbitcnt,   :mpz_remove,               [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_fac_ui,               [:M,     :UI                    ]),
    (:Cvoid,     :mpz_2fac_ui,              [:M,     :UI                    ]),
    (:Cvoid,     :mpz_mfac_uiui,            [:M,     :UI,    :UI            ]),
    (:Cvoid,     :mpz_primorial_ui,         [:M,     :UI                    ]),
    (:Cvoid,     :mpz_bin_ui,               [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_bin_uiui,             [:M,     :UI,    :UI            ]),
    (:Cvoid,     :mpz_fib_ui,               [:M,     :UI                    ]),
    (:Cvoid,     :mpz_fib2_ui,              [:M,     :M,     :UI            ]),
    (:Cvoid,     :mpz_lucnum_ui,            [:M,     :UI                    ]),
    (:Cvoid,     :mpz_lucnum2_ui,           [:M,     :M,     :UI            ]),

    # 5.10 Comparison Functions
    (:Cint,      :mpz_cmp,                  [:M,     :M                     ]),
    (:Cint,      :mpz_cmp_d,                [:M,     :D                     ]),
    (:Cint,      :mpz_cmp_si,               [:M,     :SI                    ]),
    (:Cint,      :mpz_cmp_ui,               [:M,     :UI                    ]),
    (:Cint,      :mpz_cmpabs,               [:M,     :M                     ]),
    (:Cint,      :mpz_cmpabs_d,             [:M,     :D                     ]),
    (:Cint,      :mpz_cmpabs_ui,            [:M,     :UI                    ]),

    # 5.11 Logical and Bit Manipulation Functions
    (:Cvoid,     :mpz_and,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_ior,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_xor,                  [:M,     :M,     :M             ]),
    (:Cvoid,     :mpz_com,                  [:M,     :M                     ]),
    (:Cbitcnt,   :mpz_popcount,             [:M                            ]),
    (:Cbitcnt,   :mpz_hamdist,              [:M,     :M                     ]),
    (:Cbitcnt,   :mpz_scan0,                [:M,     :BC                    ]),
    (:Cbitcnt,   :mpz_scan1,                [:M,     :BC                    ]),
    (:Cvoid,     :mpz_setbit,               [:M,     :BC                    ]),
    (:Cvoid,     :mpz_clrbit,               [:M,     :BC                    ]),
    (:Cvoid,     :mpz_combit,               [:M,     :BC                    ]),
    (:Cint,      :mpz_tstbit,               [:M,     :BC                    ]),

    # 5.12 Input and Output Functions. FILE* is intentionally represented as Ptr{Cvoid}.
    (:Csize_t,   :mpz_out_str,              [:PV,    :I,     :M             ]),
    (:Csize_t,   :mpz_inp_str,              [:M,     :PV,    :I             ]),
    (:Csize_t,   :mpz_out_raw,              [:PV,    :M                     ]),
    (:Csize_t,   :mpz_inp_raw,              [:M,     :PV                    ]),

    # 5.13 Random Number Functions. gmp_randstate_t is represented as Ptr{Cvoid}.
    (:Cvoid,     :mpz_urandomb,             [:M,     :RAND,  :BC            ]),
    (:Cvoid,     :mpz_urandomm,             [:M,     :RAND,  :M             ]),
    (:Cvoid,     :mpz_rrandomb,             [:M,     :RAND,  :BC            ]),
    (:Cvoid,     :mpz_random,               [:M,     :MS                    ]),
    (:Cvoid,     :mpz_random2,              [:M,     :MS                    ]),

    # 5.14 Integer Import and Export
    (:Cvoid,     :mpz_import,               [:M,     :SZ,    :I,     :SZ,    :I,     :SZ,    :PV]),
    (:PV,        :mpz_export,               [:PV,    :PSZ,   :I,     :SZ,    :I,     :SZ,    :M ]),

    # 5.15 Miscellaneous Functions
    (:Cint,      :mpz_fits_ulong_p,         [:M                            ]),
    (:Cint,      :mpz_fits_slong_p,         [:M                            ]),
    (:Cint,      :mpz_fits_uint_p,          [:M                            ]),
    (:Cint,      :mpz_fits_sint_p,          [:M                            ]),
    (:Cint,      :mpz_fits_ushort_p,        [:M                            ]),
    (:Cint,      :mpz_fits_sshort_p,        [:M                            ]),
    (:Csize_t,   :mpz_sizeinbase,           [:M,     :I                     ]),

    # 5.16 Special Functions
    (:Csize_t,   :mpz_size,                 [:M                            ]),
    (:mp_limb_t, :mpz_getlimbn,             [:M,     :MS                    ]),
    (:PML,       :mpz_limbs_read,           [:M                            ]),
    (:PML,       :mpz_limbs_write,          [:M,     :MS                    ]),
    (:PML,       :mpz_limbs_modify,         [:M,     :MS                    ]),
    (:Cvoid,     :mpz_limbs_finish,         [:M,     :MS                    ]),
]

function gmp_symbol(fname::Symbol)
    Symbol(:__g, fname)
end

function has_gmp_symbol(fname::Symbol)
    ptr = dlsym_e(dlopen(libgmp), gmp_symbol(fname))
    return ptr != C_NULL
end

gmp_arg(::Val, x) = x
gmp_arg(::Val{:UI}, x) = Culong(x)
gmp_arg(::Val{:SI}, x) = Clong(x)
gmp_arg(::Val{:I}, x) = Cint(x)
gmp_arg(::Val{:BC}, x) = bitcnt_t(x)
gmp_arg(::Val{:SZ}, x) = Csize_t(x)
gmp_arg(::Val{:ML}, x) = mp_limb_t(x)
gmp_arg(::Val{:MS}, x) = mp_size_t(x)

function gmp_ffi_defining_functions()
    for (ret_type, fname, arg_keys) in GMP_DEFINITIONS
        if !has_gmp_symbol(fname)
            @warn "GMP symbol $(gmp_symbol(fname)) not found. Skipping auto-definition."
            continue
        end

        j_types = [TYPE_MAP[k][1] for k in arg_keys]
        c_types = [TYPE_MAP[k][2] for k in arg_keys]
        args = [Symbol(:a, i) for i in eachindex(arg_keys)]
        c_args = [:($(gmp_arg)(Val{$(QuoteNode(arg_keys[i]))}(), $(args[i]))) for i in eachindex(args)]
        sig = [Expr(:(::), args[i], j_types[i]) for i in eachindex(args)]
        j_ret = RET_MAP[ret_type][1]
        c_ret = RET_MAP[ret_type][2]
        target = (gmp_symbol(fname), libgmp)

        @eval begin
            export $fname
            function $fname($(sig...))::$j_ret
                return ccall($target, $c_ret, ($(c_types...),), $(c_args...))
            end
        end
    end
end

gmp_ffi_defining_functions()

export mpz_mod_ui, mpz_kronecker, mpz_sgn, mpz_odd_p, mpz_even_p

# GMP exposes these as macros or compatibility aliases, not as exported DLL symbols.
mpz_mod_ui(a1::BigInt, a2::BigInt, a3::CulongArg)::Culong = mpz_fdiv_r_ui(a1, a2, a3)
mpz_kronecker(a1::BigInt, a2::BigInt)::Cint = mpz_jacobi(a1, a2)
mpz_sgn(a1::BigInt)::Cint = Cint(sign(mpz_cmp_ui(a1, Culong(0))))
mpz_odd_p(a1::BigInt)::Cint = mpz_tstbit(a1, Culong(0))
mpz_even_p(a1::BigInt)::Cint = Cint(mpz_tstbit(a1, Culong(0)) == 0)
