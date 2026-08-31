import xtrack as xt
import xpart as xp
import numpy as np

'''
Whole-element CLIC Pre-Damping Ring (PDR)

The PDR has a racetrack geometry, with TME cells in the arcs and FODOW
cells (FODO cells containing wigglers) in the straight sections. The arcs
and straights are connected through two distinct dispersion-suppression
and beta-matching sections: a DS section from arc to straight and an MDS
section from straight to arc.

Two wiggler models are available: a hard-edge model constructed from
xt.RBend wiggler pole elements and drifts and a sinusoidal planar wiggler field model implemented using
xtrack.SplineBoris.

Whole elements are retained to keep element shifts and fringe-field
geometry consistent with the physical PDR. Symmetric cell matching (fodow's and tme's) is performed
using shallow copies of these elements, with each matching section beginning
and ending at a symmetry point, yielding the appropriate lattice-design knobs.
'''


#########################################
# clic pdr beam definition #
#########################################
particle = "positron"               
energy0_gev = 2.86                                      
npart_per_bunch = 4.5e9   ## n/a          
kbunch = 156   ## n/a 
n_trains = 2   ## n/a 
h = 1298 # frf is 1 ghz 
ex_inj = 1.250699e-6   ## n/a                
ey_inj = 1.250699e-3   ## n/a               
et = 7.5e-5   ## n/a                      
sigt_m = 0.009   ## n/a 
sige = 0.015   ## n/a 

particle_ref = xp.Particles(
    mass0=xp.ELECTRON_MASS_EV,        
    q0=+1,                            
    energy0=energy0_gev * 1e9          
)


#########################################
# xsuite clic pdr #
#########################################   
env = xt.get_environment()            
env.particle_ref=particle_ref
env.vars.default_to_zero = False              ###

#########################################
# rigidity
#########################################
env.vars.update({  
    "nd": 38,           # number of dipoles [-]
    'ntme': 32,         # number of tme cells
    "nfodow": 18,       # number of fodow cells
    "en": 2.86,         # energy [GeV]
    "bb": 1.2,          # arc bending field [T]
    "bw": 1.9,          # wiggler field [T]  ## for the hard-edge wiggler model
})

env['trep'] = 20 # repetition rate [ms] 

#########################################
# tme parameters
#########################################

# tme drift lengths
env.vars.update({
    "l1a": 0.9,
    "l2a": 0.6,
    "l3a": 0.5 * 2,
    "lls": 0.3,  # length of all sexts
})

#########################################
# arc cell: tme elements
#########################################

# ---dipole--- 
env['nd'] = 38
env["ang"] = 2*np.pi / (env["nd"])     # bending angle of whole dipole [rad]
env["r"]   = (1/0.2998) * (env["en"] / env["bb"])
env["s"]  = env["ang"] * env["r"]                 # dipole length [m] = 1.31
env["ld"] = "ang * r"

env.new(
    "dip", xt.Bend,
    length=env["ld"], angle=env["ang"],
)

# ---arc quads---
env.vars.update({  
    "kqfa":  2.492391737,  ### initial values for the matching
    "kqda": -2.069416029,
    "lq": 0.28,
})

env.new("qfa", xt.Quadrupole, length=env["lq"],   k1="kqfa")
env.new("qda", xt.Quadrupole, length=env["lq"],   k1="kqda")

# ---arc sexts---
env.vars.update({
    "ksx1": 0.0,
    "ksy1": 0.0,
})

env.new("sx1", xt.Sextupole,  length=env["lls"], k2="ksx1")
env.new("sy1", xt.Sextupole,  length=env["lls"], k2="ksy1")

#########################################
# tme cell
#########################################

##### define the edge to edge drift lengths #####
env.vars.update({
    "dtme_1": (env["l1a"] - env["lls"]) / 2,
    "dtme_2": (env["l1a"] - env["lls"] - env["lq"]) / 2,
    "dtme_3": env["l2a"] - env["lq"],
    "dtme_4": (env["l3a"] - env["lq"] - env["lls"]) / 2,
})

# drift elements
env.new("d_tme_1", xt.Drift, length=env["dtme_1"]) # dip2->sx1
env.new("d_tme_2", xt.Drift, length=env["dtme_2"]) # sx1->qfa
env.new("d_tme_3", xt.Drift, length=env["dtme_3"]) # qfa->qda
env.new("d_tme_4", xt.Drift, length=env["dtme_4"]) # qda->sy1

# ------tme build------
env.new_line(
    name="tme",
    components=[
    'dip',
    env.new('d_tme_1.1', 'd_tme_1'),
    env.new('sx1.1', 'sx1'),        ## element.copy() brings the knobs
    env.new('d_tme_2.1', 'd_tme_2'),
    env.new('qfa.1', 'qfa'),
    env.new('d_tme_3.1', 'd_tme_3'),
    env.new('qda.1', 'qda'),
    env.new('d_tme_4.1', 'd_tme_4'),
    env.new('sy1.1', 'sy1'), ### symmetry point
    env.new('d_tme_4.2', 'd_tme_4'),
    env.new('qda.2', 'qda'),
    env.new('d_tme_3.2', 'd_tme_3'),
    env.new('qfa.2', 'qfa'),
    env.new('d_tme_2.2', 'd_tme_2'),
    env.new('sx1.2', 'sx1'),
    env.new('d_tme_1.2', 'd_tme_1'),\
],
)

#########################################
# wiggler parameters
#########################################

# wiggler lengths and bending angle
env['nwigg'] = 36
env["wigglerperiod"] = 0.3                       # [m]
env["polelength"] = env["wigglerperiod"] / 4.0   # [m]
env["wigangle"] = env["polelength"] * env["bw"] / (3.356 * env["en"])  # wiggler bend [rad]

#########################################
# wiggler poles (hard-edge model)
#########################################

## Rbend defined by k0
env["wigg_k0_pos"] = env["wigangle"] / env["polelength"]
env["wigg_k0_neg"] = -env["wigangle"] / env["polelength"]

env.vars.update({ 
    "wigg_k0_pos": env["wigg_k0_pos"],
    "wigg_k0_neg": env["wigg_k0_neg"],
})

env["wigg_fdown_pos"] =  env["wigangle"] / 2.0 # angle of the reference trajectory at the edge: to produce rectangular y-focusing
env["wigg_fdown_neg"] = -env["wigangle"] / 2.0

env["wigg_fdown_halfpos"] =  env["wigangle"] / 4.0
env["wigg_fdown_halfneg"] = -env["wigangle"] / 4.0

# RBend defined by k0
### This has the reference trajectory follow the racetrack geometry, and the closed orbit wiggle around it.
wigpolepos = env.new("wigpolepos", xt.RBend,length_straight=env["polelength"],
    angle=0.0,
    k0="wigg_k0_pos",
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_angle_fdown="wigg_fdown_pos",  ##
    edge_exit_angle_fdown="wigg_fdown_pos",   ##
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

wigpoleneg = env.new(
    "wigpoleneg", xt.RBend,
    length_straight=env["polelength"],
    angle=0.0,
    k0="wigg_k0_neg",
    k0_from_h=False,
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_angle_fdown="wigg_fdown_neg",
    edge_exit_angle_fdown="wigg_fdown_neg",
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

wighalfpolepos = env.new(
    "wighalfpolepos", xt.RBend,
    length_straight=env["polelength"] / 2.0,
    angle=0.0,
    k0="wigg_k0_pos",
    k0_from_h=False,
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_angle_fdown="wigg_fdown_halfpos",
    edge_exit_angle_fdown="wigg_fdown_halfpos",
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

wighalfpoleneg = env.new(
    "wighalfpoleneg", xt.RBend,
    length_straight=env["polelength"] / 2.0,
    angle=0.0,
    k0="wigg_k0_neg",
    k0_from_h=False,
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_angle_fdown="wigg_fdown_halfneg",
    edge_exit_angle_fdown="wigg_fdown_halfneg",
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

'''
# RBend defined by angle (madx way)
### This bends the reference trajectory itself to follow the particle path
wigpolepos = env.new(
    "wigpolepos", xt.RBend,
    length_straight=env["polelength"],
    angle=env["wigangle"],
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

wigpoleneg = env.new(
    "wigpoleneg", xt.RBend,
    length_straight=env["polelength"],
    angle=-env["wigangle"],
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

wighalfpolepos = env.new(
    "wighalfpolepos", xt.RBend,
    length_straight=env["polelength"] / 2.0,
    angle=env["wigangle"] / 2.0,
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

wighalfpoleneg = env.new(
    "wighalfpoleneg", xt.RBend,
    length_straight=env["polelength"] / 2.0,
    angle=-env["wigangle"] / 2.0,
    edge_entry_angle=0.0,
    edge_exit_angle=0.0,
    edge_entry_model="linear",
    edge_exit_model="linear",
    rbend_model="straight-body",
    force=True,
)

'''

# drift between poles  
wigdrift = env.new("wigdrift", xt.Drift, length=env["polelength"])

#########################################
# wiggler build
#########################################

wig_components = [wighalfpolepos, wigdrift]  # START with one halfpolepos and drift 

for i in range(19):           # -> then alternating between 19 poles (10 poleneg and 9 polepos) including drifts inbetween
    wig_components.append(wigpoleneg if (i % 2 == 0) else wigpolepos)
    wig_components.append(wigdrift)

wig_components.append(wighalfpolepos)  # -> FINISH by adding one halfpolepos

env.new_line(name="wiggler", components=wig_components)
env["wiggler"].replace_all_repeated_elements() ###

#########################################
######### wiggler sinusoidal ############
#########################################

## modelling the sinusoidal planar wiggler with xtrack.SplineBoris ##

#########################################
# 1. wiggler design parameters
#########################################

Bw = 1.9   # peak vertical field [T]
lambda_w = 0.30   # wiggler period [m]
n_periods = 10    # total number of periods (10 pos and 10 neg poles)
k = 2*np.pi/lambda_w

L_wig = n_periods * lambda_w

splines_per_period = 8 ### spline fit choice: one spline for each half-pole
n_splines = n_periods * splines_per_period
L_spline = L_wig / n_splines

# Boris integration steps per spline piece.
# This affects tracking accuracy, not the spline fit itself.
n_splineboris = 10

##############################################################
# 2. analytic on-axis BPMeth coefficient
##############################################################

def b1(s):
    return Bw * np.cos(k*s)

def db1_ds(s):
    return -Bw * k * np.sin(k*s)

##############################################################
# 3. build one Spline4 for an interval ds
##############################################################

from scipy.integrate import quad

def make_spline4_from_function(f, df, s0, s1):
    L = s1 - s0
    integ_avg = quad(f, s0, s1, epsabs=1e-14, epsrel=1e-14)[0] / L

    return xt.Spline4(
        val_start=f(s0),
        der_start=df(s0),
        val_end=f(s1),
        der_end=df(s1),
        mean=integ_avg,
    )

def zero_spline4(s0, s1):
    return xt.Spline4(
        val_start=0.0,
        der_start=0.0,
        val_end=0.0,
        der_end=0.0,
        mean=0.0,
    )

##############################################################
# 4. splineBoris wigglers
##############################################################

def make_sb_wiggler_line(line_name):
    sb_components = []

    for ii in range(n_splines):
        s0 = ii * L_spline
        s1 = (ii + 1) * L_spline

        by0 = make_spline4_from_function(b1, db1_ds, s0, s1)
        bs0 = zero_spline4(s0, s1)

        name = f"{line_name}.sb_{ii:03d}"

        env.elements[name] = xt.SplineBoris(
            bs=bs0,
            by=(by0,),
            bx=(None,),
            length=L_spline,
            n_steps=n_splineboris,
        )

        sb_components.append(name)

    env.new_line(name=line_name, components=sb_components)

    return env[line_name]

#########################################
# f0d0w elements
#########################################

# f0d0w lengths
env["lq1w"] = 0.1 * 2      # [m] quad length
env["lq2w"] = 0.1 * 2      # [m]
env["lwd1"] = 0.5          # [m] drift length
env["lwd2"] = 0.5          # [m]

env.vars.update({               
    "kqfw":  1.55261226,
    "kqdw": -1.301643467,
})

env.new("qfw", xt.Quadrupole, length=env["lq1w"], k1="kqfw")  
env.new("qdw", xt.Quadrupole, length=env["lq2w"], k1="kqdw")

env.new("wd1",  xt.Drift, length=env["lwd1"])
env.new("wd2",  xt.Drift, length=env["lwd2"])

# -----fodow build-----
# Hard-edge wiggler fodow:
env.new_line(
    name="fodow",
    components=[
        'qfw',
        env.new('wd2.1', 'wd2'),
        env["wiggler"].clone(suffix=".1"),
        env.new('wd1.1', 'wd1'),
        'qdw',
        env.new('wd1.2', 'wd1'),
        env["wiggler"].clone(suffix=".2"),
        env.new('wd2.2', 'wd2'),
    ],
)

# SplineBoris wigglers fodow:
make_sb_wiggler_line("wiggler.1_sb")
make_sb_wiggler_line("wiggler.2_sb")

env.new_line(
    name="fodow_sb",
    components=[
        env.new('qfw_sb', 'qfw'),
        env.new('wd2.1_sb', 'wd2'),
        env["wiggler.1_sb"],
        env.new('wd1.1_sb', 'wd1'),
        env.new('qdw_sb', 'qdw'),
        env.new('wd1.2_sb', 'wd1'),
        env["wiggler.2_sb"],
        env.new('wd2.2_sb', 'wd2'),
    ],
)

#########################################
# dispsup-beta matching sections
#########################################

# ds & mds quads
env["lqds"]   = 0.35

env.vars.update({
    "kqds11": -0.951299623,     ### ds knobs: kqds11-17         
    "kqds12":  3.215394621,
    "kqds13": -2.322720333,
    "kqds14": -1.22300104,
    "kqds15":  2.078486212,
    "kqds16": -0.1274254316,
    "kqds17": -1.306390283,     
    "kqds18":  0.8174937228,    ### mds knobs: kqds18,19 + kqds21-17
    "kqds19": -0.7202944284,
    "kqds21":  2.19351565,
    "kqds22": -3.193124933,
    "kqds23":  2.952413639,
    "kqds24": -1.518275912,
    "kqds25":  0.5014376022,
    "kqds26":  1.323682568,
    "kqds27": -1.411675205,
}) 

# ds quads
env.new("qds11",  xt.Quadrupole, length=env["lqds"], k1="kqds11")
env.new("qds12",  xt.Quadrupole, length=env["lqds"], k1="kqds12")
env.new("qds13",  xt.Quadrupole, length=env["lqds"], k1="kqds13")
env.new("qds14",  xt.Quadrupole, length=env["lqds"], k1="kqds14")
env.new("qds15",  xt.Quadrupole, length=env["lqds"], k1="kqds15")
env.new("qds16",  xt.Quadrupole, length=env["lqds"], k1="kqds16")
env.new("qds17",  xt.Quadrupole, length=env["lqds"], k1="kqds17")
env.new("qds18",  xt.Quadrupole, length=env["lqds"], k1="kqds18")
env.new("qds19",  xt.Quadrupole, length=env["lqds"], k1="kqds19")

# mds quads
env.new("qds21", xt.Quadrupole, length=env["lqds"], k1="kqds21")
env.new("qds22", xt.Quadrupole, length=env["lqds"], k1="kqds22")
env.new("qds23", xt.Quadrupole, length=env["lqds"], k1="kqds23")
env.new("qds24", xt.Quadrupole, length=env["lqds"], k1="kqds24")
env.new("qds25", xt.Quadrupole, length=env["lqds"], k1="kqds25")
env.new("qds26", xt.Quadrupole, length=env["lqds"], k1="kqds26")
env.new("qds27", xt.Quadrupole, length=env["lqds"], k1="kqds27")

# ds & mds sext strengths (N/A)
env.vars.update({
    "ksx2": 0.0,     
    "ksy2": 0.0,
})

env.new("sx2", xt.Sextupole,  length=env["lls"], k2="ksx2")
env.new("sy2", xt.Sextupole,  length=env["lls"], k2="ksy2")

# all ds-bm regions (ds & mds) drift lengths
env.vars.update({
    "lds11": 0.6,
    "lds12": 0.6,
    "lds13": 0.5,
    "lds14": 0.5,
    "lds15": 0.5,
    "lds16": 0.5,
    "lds17": 2.0,
    "lds18": 2.0,# mds_2
    "lds19": 2.0,# mds_1
    "lds110": 3.0,# ds_10
    "lds111": 3.0,# ds_11

    "lds21": 0.6,# ds_1 # mds_9
    "lds22": 0.6,# ds_2 # mds_8
    "lds23": 0.5,# ds_3 # mds_7
    "lds24": 0.5,# ds_4 # mds_6
    "lds25": 0.4,# ds_5 # mds_5
    "lds26": 0.6,# ds_6 # mds_4
    "lds27": 0.6,# ds_7 # mds_3
    "lds28": 3.3,# ds_8  
    "lds29": 3.3 ,# ds_9 
})

# ds total length
env["lds"] = (
    env["ld"]
  + env["lds21"] + env["lqds"]
  + env["lds22"] + env["lqds"]
  + env["lds23"] + env["lqds"]
  + env["lds24"]
  + env["ld"]
  + env["lds25"] + env["lqds"]
  + env["lds26"] + env["lqds"]
  + env["lds27"] + env["lqds"]
  + env["lds28"] + env["lqds"]
  + env["lds29"] + env["lqds"]
  + env["lds110"] + env["lqds"]
  + env["lds111"]
)

# ----------ds build----------
# define end Drift
env.new("dend", xt.Drift, length=env["lds111"])

# ds_1
env.new_line(
    name="ds_1",
    length=env["lds"],
    components=[
        env.place(env.new("dip.ds_1.1", "dip"), anchor="center", at=env["ld"]/2),

        env.place(env.new("sx2.ds_1.1", "sx2"), anchor="center",
                  at=env["ld"] + (env["lds21"]-env["lls"])/2 + env["lls"]/2),

        env.place(env.new("qds11.ds_1.1", "qds11"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]/2),

        env.place(env.new("qds12.ds_1.1", "qds12"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]/2),

        env.place(env.new("sy2.ds_1.1", "sy2"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + (env["lds23"]-env["lls"])/2 + env["lls"]/2),

        env.place(env.new("qds13.ds_1.1", "qds13"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]/2),

        env.place(env.new("dip.ds_1.2", "dip"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]/2),

        env.place(env.new("qds14.ds_1.1", "qds14"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]/2),

        env.place(env.new("qds15.ds_1.1", "qds15"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]/2),

        env.place(env.new("qds16.ds_1.1", "qds16"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]/2),

        env.place(env.new("qds17.ds_1.1", "qds17"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]/2),

        env.place(env.new("qds18.ds_1.1", "qds18"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]
                  + env["lds29"] + env["lqds"]/2),

        env.place(env.new("qds19.ds_1.1", "qds19"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]
                  + env["lds29"] + env["lqds"]
                  + env["lds110"] + env["lqds"]/2),

        env.place(env.new("dend.ds_1.1", "dend"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]
                  + env["lds29"] + env["lqds"]
                  + env["lds110"] + env["lqds"]
                  + env["lds111"]/2),
    ],
)

# ds_2
env.new_line(
    name="ds_2",
    length=env["lds"],
    components=[
        env.place(env.new("dip.ds_2.1", "dip"), anchor="center", at=env["ld"]/2),

        env.place(env.new("sx2.ds_2.1", "sx2"), anchor="center",
                  at=env["ld"] + (env["lds21"]-env["lls"])/2 + env["lls"]/2),

        env.place(env.new("qds11.ds_2.1", "qds11"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]/2),

        env.place(env.new("qds12.ds_2.1", "qds12"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]/2),

        env.place(env.new("sy2.ds_2.1", "sy2"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + (env["lds23"]-env["lls"])/2 + env["lls"]/2),

        env.place(env.new("qds13.ds_2.1", "qds13"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]/2),

        env.place(env.new("dip.ds_2.2", "dip"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]/2),

        env.place(env.new("qds14.ds_2.1", "qds14"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]/2),

        env.place(env.new("qds15.ds_2.1", "qds15"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]/2),

        env.place(env.new("qds16.ds_2.1", "qds16"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]/2),

        env.place(env.new("qds17.ds_2.1", "qds17"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]/2),

        env.place(env.new("qds18.ds_2.1", "qds18"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]
                  + env["lds29"] + env["lqds"]/2),

        env.place(env.new("qds19.ds_2.1", "qds19"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]
                  + env["lds29"] + env["lqds"]
                  + env["lds110"] + env["lqds"]/2),

        env.place(env.new("dend.ds_2.1", "dend"), anchor="center",
                  at=env["ld"] + env["lds21"] + env["lqds"]
                  + env["lds22"] + env["lqds"]
                  + env["lds23"] + env["lqds"]
                  + env["lds24"] + env["ld"]
                  + env["lds25"] + env["lqds"]
                  + env["lds26"] + env["lqds"]
                  + env["lds27"] + env["lqds"]
                  + env["lds28"] + env["lqds"]
                  + env["lds29"] + env["lqds"]
                  + env["lds110"] + env["lqds"]
                  + env["lds111"]/2),
    ],
)
# rf cavity in ds
env["lrf"]  = 1e-5 # [m] ### thin cavity
env['h']=2596/2
env.new("rf", xt.Cavity, length=env["lrf"], voltage=10e6, frequency=0.0, lag=0, force=True, harmonic=env['h'])  # 10 mev

# mds total length
env["lds2"] = (
    env["lq1w"]
  + env["lds19"] + env["lqds"]
  + env["lds18"] + env["lqds"]
  + env["lds27"] + env["lqds"]
  + env["lds26"] + env["lqds"]
  + env["lds25"]
  + env["ld"]
  + env["lds24"] + env["lqds"]
  + env["lds23"] + env["lqds"]
  + env["lds22"] + env["lqds"]
  + env["lds21"]
)

# ----------mds build----------
# define end drift
env.new("dend2", xt.Drift, length=(env["lds21"]-env['lls'])/2)

# mds_1
env.new_line(
    name="mds_1",
    length=env["lds2"],
    components=[

        env.place(env.new("qfw.mds_1.1", "qfw"),
            anchor="center",
            at=env["lq1w"]/2),

        env.place(env.new("qds27.mds_1.1", "qds27"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]/2),

        env.place(env.new("qds26.mds_1.1", "qds26"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]/2),

        env.place(env.new("qds25.mds_1.1", "qds25"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]/2),

        env.place(env.new("qds24.mds_1.1", "qds24"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]/2),

        env.place(env.new("dip.mds_1.1", "dip"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]/2),

        env.place(env.new("qds23.mds_1.1", "qds23"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]/2),

        env.place(env.new("sy2.mds_1.1", "sy2"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + (env["lds23"] - env["lls"])/2 + env["lls"]/2),

        env.place(env.new("qds22.mds_1.1", "qds22"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]/2),

        env.place(env.new("qds21.mds_1.1", "qds21"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]
              + env["lds22"] + env["lqds"]/2),

        env.place(env.new("sx2.mds_1.1", "sx2"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]
              + env["lds22"] + env["lqds"]
              + (env["lds21"] - env["lls"])/2 + env["lls"]/2),

        env.place(env.new("dend2.mds_1.1", "dend2"),
            anchor="start",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]
              + env["lds22"] + env["lqds"]
              + (env["lds21"] + env["lls"])/2),
    ],
)

# mds_2
env.new_line(
    name="mds_2",
    length=env["lds2"],
    components=[

        env.place(env.new("qfw.mds_2.1", "qfw"),
            anchor="center",
            at=env["lq1w"]/2),

        env.place(env.new("qds27.mds_2.1", "qds27"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]/2),

        env.place(env.new("qds26.mds_2.1", "qds26"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]/2),

        env.place(env.new("qds25.mds_2.1", "qds25"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]/2),

        env.place(env.new("qds24.mds_2.1", "qds24"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]/2),

        env.place(env.new("dip.mds_2.1", "dip"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]/2),

        env.place(env.new("qds23.mds_2.1", "qds23"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]/2),

        env.place(env.new("sy2.mds_2.1", "sy2"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + (env["lds23"] - env["lls"])/2 + env["lls"]/2),

        env.place(env.new("qds22.mds_2.1", "qds22"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]/2),

        env.place(env.new("qds21.mds_2.1", "qds21"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]
              + env["lds22"] + env["lqds"]/2),

        env.place(env.new("sx2.mds_2.1", "sx2"),
            anchor="center",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]
              + env["lds22"] + env["lqds"]
              + (env["lds21"] - env["lls"])/2 + env["lls"]/2),

        env.place(env.new("dend2.mds_2.1", "dend2"),
            anchor="start",
            at=env["lq1w"]
              + env["lds19"] + env["lqds"]
              + env["lds18"] + env["lqds"]
              + env["lds27"] + env["lqds"]
              + env["lds26"] + env["lqds"]
              + env["lds25"] + env["ld"]
              + env["lds24"] + env["lqds"]
              + env["lds23"] + env["lqds"]
              + env["lds22"] + env["lqds"]
              + (env["lds21"] + env["lls"])/2),
    ],
)

############ build pdr sections with whole element cells 
#  & select symmetric cells with markers for the matching ############

 # ARCS
cells = []
for ii in range(env["ntme"]):     # order all TME cells inside a circular ring
    new_cell = env['tme'].clone(suffix=f'tme_{ii+1}')
    cells.append(new_cell)

env.new_line(
    name="ring",
    components=cells)

## prepare for matching: insert markers to make a symmetric cell ##
ring = env['ring'].copy(shallow=True)

tt = ring.get_table()     # create copy of ring to put markers in

tt_bend = tt.rows.match(element_type='Bend')

insertions = []        # insert markers at dip centers for matching
for nn in tt_bend.name:
    nn_marker = nn + '.center'
    env.new(nn_marker, xt.Marker)
    insertions.append(env.place(nn_marker, at=0, from_=nn)) #, from_anchor='center'))

# symmetric tme cell
ring.insert(insertions)
env["tme_cell_s"] = ring.select('dip.tme_8.center', 'dip.tme_9.center') # match this: "tme_cell_s"

env["arc1"] = ring.select(    # divide circular ring into 2 arcs (shallow copies with whole elements)
    start="dip.tme_1_entry",  # these lines have the same knobs as the "ring" line. Check with: env.info("kqfa", limit=None)
    end="dip.tme_17_entry",
)

env["arc2"] = ring.select(
    start="dip.tme_17_entry",
    end='_end_point',
)

# STRAIGHTS
# hard-edge wiggler straight:
cells = []
for ii in range(env["nfodow"]):
    new_cell = env['fodow'].clone(suffix=f'fodow_{ii+1}')
    cells.append(new_cell)

env.new_line(
    name="straight",
    components=cells)

straight = env['straight'].copy(shallow=True)

tt = straight.get_table()

tt_quad = tt.rows.match(element_type='Quadrupole')

insertions = []
for nn in tt_quad.name:
    nn_marker = nn + '.center'
    env.new(nn_marker, xt.Marker)
    insertions.append(env.place(nn_marker, at=0, from_=nn)) #, from_anchor='center'))

# symmetric hard-edge fodow cell
straight.insert(insertions)
env["fodow_cell_s"] = straight.select('qfw.fodow_3.center', 'qfw.fodow_4.center') # match this symmetric cell

env["ss1"] = straight.select(
    start="qfw.fodow_1_entry",
    end="qfw.fodow_10_entry",
)

env["ss2"] = straight.select(
    start="qfw.fodow_10_entry",
    end='_end_point',
)

# sinusoidal wiggler straight (SplineBoris) 
cells_sb = []
for ii in range(env["nfodow"]):
    i_cell = ii + 1

    make_sb_wiggler_line(f"wiggler.{i_cell}.1_sb")
    make_sb_wiggler_line(f"wiggler.{i_cell}.2_sb")

    env.new_line(
        name=f"fodow_{i_cell}_sb",
        components=[
            env.new(f"qfw.fodow_{i_cell}_sb", "qfw"),
            env.new(f"wd2.1.fodow_{i_cell}_sb", "wd2"),
            env[f"wiggler.{i_cell}.1_sb"],
            env.new(f"wd1.1.fodow_{i_cell}_sb", "wd1"),
            env.new(f"qdw.fodow_{i_cell}_sb", "qdw"),
            env.new(f"wd1.2.fodow_{i_cell}_sb", "wd1"),
            env[f"wiggler.{i_cell}.2_sb"],
            env.new(f"wd2.2.fodow_{i_cell}_sb", "wd2"),
        ],
    )

    cells_sb.append(env[f"fodow_{i_cell}_sb"])

env.new_line(
    name="straight_sb",
    components=cells_sb,
)

straight_sb = env['straight_sb'].copy(shallow=True)

tt_sb = straight_sb.get_table()

tt_quad_sb = tt_sb.rows.match(element_type='Quadrupole')

insertions_sb = []
for nn in tt_quad_sb.name:
    nn_marker = nn + '.center'
    env.new(nn_marker, xt.Marker)
    insertions_sb.append(env.place(nn_marker, at=0, from_=nn)) #, from_anchor='center'))

# symmetric SplineBoris fodow cell (if rematch is needed)
straight_sb.insert(insertions_sb)
env["fodow_cell_s_sb"] = straight_sb.select('qfw.fodow_3_sb.center', 'qfw.fodow_4_sb.center')

env["ss1_sb"] = straight_sb.select(
    start="qfw.fodow_1_sb_entry",
    end="qfw.fodow_10_sb_entry",
)

env["ss2_sb"] = straight_sb.select(
    start="qfw.fodow_10_sb_entry",
    end='_end_point',
)

# DS #

con_ds = env["ds_1"].copy(shallow=True) + env["ss2"].copy(shallow=True)

env.new("dip.ds_1.1.center", xt.Marker)
con_ds.insert(
    env.place("dip.ds_1.1.center", at=0, from_="dip.ds_1.1")
)

env["ds_match"] = con_ds.select(
    start="dip.ds_1.1.center",
    end="qfw.fodow_10.center",
) # match ds starting from middle of arc's dip and ending in middle of straight's quad

# MDS #

con_mds = env["mds_1"].copy(shallow=True) + env["arc1"].copy(shallow=True)

env.new("qfw.mds_1.1.center", xt.Marker)
con_mds.insert(
    env.place("qfw.mds_1.1.center", at=0, from_="qfw.mds_1.1")
)

env["mds_match"] = con_mds.select(
    start="qfw.mds_1.1.center",
    end="dip.tme_1.center",
) # match mds starting from middle of straight's quad and ending in middle of arc's dip

# ============================================================
# SplineBoris DS match line: arc -> SplineBoris wiggler straight
# ============================================================

con_ds_sb = env["ds_1"].copy(shallow=True) + env["ss2_sb"].copy(shallow=True)

env.new("dip.ds_1.1.center_sb", xt.Marker)
con_ds_sb.insert(
    env.place("dip.ds_1.1.center_sb", at=0, from_="dip.ds_1.1")
)

env["ds_match_sb"] = con_ds_sb.select(
    start="dip.ds_1.1.center_sb",
    end="qfw.fodow_10_sb.center",
)

# ============================================================
# SplineBoris MDS match line: SplineBoris straight -> MDS -> arc
# ============================================================

con_mds_sb = env["mds_1"].copy(shallow=True) + env["arc1"].copy(shallow=True)

env.new("qfw.mds_1.1.center_sb", xt.Marker)
con_mds_sb.insert(
    env.place("qfw.mds_1.1.center_sb", at=0, from_="qfw.mds_1.1")
)

env["mds_match_sb"] = con_mds_sb.select(
    start="qfw.mds_1.1.center_sb",
    end="dip.tme_1.center",
)

#########################################
# pdr build
#########################################

env["pdr"] = (  # hard-edge reference PDR
    env["ss1"]
    + env["mds_1"]
    + env["arc1"]
    + env["ds_1"]
    + env["ss2"]
    + env["mds_2"]
    + env["arc2"]
    + env["ds_2"]
)

env["pdr_sb"] = (  # SplineBoris PDR with the same non-wiggler sections/knobs
    env["ss1_sb"]
    + env["mds_1"]
    + env["arc1"]
    + env["ds_1"]
    + env["ss2_sb"]
    + env["mds_2"]
    + env["arc2"]
    + env["ds_2"]
)

env.vars.default_to_zero = False