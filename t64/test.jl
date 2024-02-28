
using LaMEM, GeophysicalModelGenerator, Plots
include("subindshapes.jl")  # Include the file containing shape definitions and functions
include("custom_gmg.jl")  # Include the file containing custom functions for GeophysicalModelGenerator

modelHeight = 2880.0km
T_mantle = 1350  # in Celcius

refDensity = 3200.0kg / m^3
deltaRhoMax = 80.0kg / m^3
gravity = 9.8m / s^2
refViscosity = 5e20Pas
bodyForce = deltaRhoMax * gravity
K_tau = bodyForce * modelHeight

CharUnits = GEO_units(length=modelHeight, stress=K_tau, viscosity=refViscosity)

modelHeight = modelHeight.val

resFactor = 2.0

model = Model(
    Grid(
        coord_x=[0, 4464, 6192, 11520],
        nel_x=round.(Int, [365, 480, 435] / resFactor),
        # nel_x=round.(Int, [365, 480*2, 435] / resFactor),

        bias_x=[0.25, 1.0, 4],
        y=[-2.5, 2.5],
        nel_y=1,
        # coord_z=[-2880, -384, 192],
        coord_z=[-2880, -384, 0],
        nel_z=round.(Int, [320, 320] / resFactor),
        bias_z=[0.25, 1.0]
    ),
    # BoundaryConditions(noslip=[0, 0, 1, 1, 0, 0]), # (left right front back bottom top)
    Scaling(CharUnits)
)

refViscosity = refViscosity.val
opViscFactor = 1.0
modelMaterials = [
    # Phase(ID=0, Name="UpperMantle", eta=refViscosity, rho=3200.0),
    Phase(ID=0, Name="UpperMantle", eta=refViscosity, eta0=refViscosity, n=3.5, e0=2.5e-15, rho=3200.0,eta_st=0.1*refViscosity,eta_vp=0.1*refViscosity),
    Phase(ID=1, Name="UppperCrustIndoAustralianPlate", eta=1e3 * refViscosity, rho=3280.0, ch=6e6), #,eta_st=0.1*refViscosity,eta_vp=0.1*refViscosity),
    Phase(ID=2, Name="LowerCrustIndoAustralianPlate", eta=1e3 * refViscosity, rho=3280.0, ch=50e6),
    Phase(ID=3, Name="LithosphericMantleCrustIndoAustralianPlate", eta=1e3 * refViscosity, rho=3280.0, ch=350e6),
    Phase(ID=4, Name="LithosphericMantleIndoAustralianPlate", eta=5e1 * refViscosity, rho=3280.0, ch=30e6),
    Phase(ID=5, Name="UpperCrustIndianIndentor", eta=1e3 * refViscosity, rho=2800.0, ch=06e6),
    Phase(ID=6, Name="UpperCrustIndianIndentor", eta=1e3 * refViscosity, rho=2800.0, ch=100e6),
    Phase(ID=7, Name="WeakULCrustIndianIndentor", eta=1e3 * refViscosity, rho=3000.0, ch=06e6),
    Phase(ID=8, Name="LowerCrustIndianIndentor", eta=1e3 * refViscosity, rho=3100.0, ch=200e6),
    Phase(ID=9, Name="WeakLLCrustIndianIndentor", eta=1e3 * refViscosity, rho=3100.0, ch=06e6),
    Phase(ID=10, Name="LithosphericMantleIndianIndentor", eta=1e3 * refViscosity, rho=3100.0, ch=350e6),
    Phase(ID=11, Name="LithosphericMantleIndianIndentor", eta=5e1 * refViscosity, rho=3200.0, ch=30e6),
    Phase(ID=12, Name="CrustEurasianPlateForeArc", eta=1e3 * refViscosity * opViscFactor, rho=2800.0),
    Phase(ID=13, Name="LithosphericMantleEurasianPlateForeArc", eta=5e2 * refViscosity * opViscFactor, rho=3220.0),
    Phase(ID=14, Name="CrustEurasianPlateBackArc", eta=2.5e2 * refViscosity * opViscFactor, rho=2800.0),
    Phase(ID=15, Name="LithosphericMantleEurasianPlateBackArc", eta=1.25e2 * refViscosity * opViscFactor, rho=3220.0),
    Phase(ID=16, Name="CrustEurasianPlate", eta=2e3 * refViscosity * opViscFactor, rho=2800.0),
    Phase(ID=17, Name="LithosphericMantleEurasianPlate", eta=1e3 * refViscosity * opViscFactor, rho=3220.0),
]


slabshapes = Slab2D(
    top_x=2082.0, top_y=0.0, length=3672.0, taper=30.0, dip=29.0, depth=180.0,
    thickness_array=[10.0, 15.0, 30.0, 30.0])

indentorshapes = Indentor2D(
    top_x=864.0, top_y=0.0, length=2448.0, taper=18.0, taper2=12.0,
    thickness_array=[5.0, 7.5, 7.5, 12.5, 7.5, 30.0, 30.0,])

overRidingShapesForeArc = OverRidingPlate2D(
    top_x=11520.0, top_y=0.0, length=-5760.0, taper=90.0, dip=29.0,
    thickness_array=[30.0, 30.0])

overRidingShapesBackArc = OverRidingPlate2D(
    top_x=11520.0, top_y=0.0, length=-5586.0, taper=90.0, dip=90.0,
    thickness_array=[30.0, 30.0])

overRidingShapes = OverRidingPlate2D(
    top_x=11520.0, top_y=0.0, length=-5040.0, taper=90.0, dip=17.0,
    thickness_array=[40.0, 80.0])

modelshapePolygons = vcat(
    slabshapes.polygons,
    indentorshapes.polygons,
    overRidingShapesForeArc.polygons,
    overRidingShapesBackArc.polygons,
    overRidingShapes.polygons)

add_phase!(model, modelMaterials[1])
model.Grid.Phases .= 0;

for (index, polygon) in enumerate(modelshapePolygons)
    AddPoly!(model, vertices=polygon, phase=ConstantPhase(index))
    add_phase!(model, modelMaterials[index+1])
end



model.Solver = Solver(
    SolverType="direct",
    DirectSolver="mumps",
    DirectPenalty=1e7,
    # SolverType      = "multigrid",
    # MGLevels        = 4,
    # MGCoarseSolver 	= "mumps",
    # DirectSolver 	= "mumps",
    PETSc_options=[
        "-snes_ksp_ew",
        "-snes_ksp_ew_rtolmax 1e-4",
        "-snes_rtol 1e-3",
        "-snes_atol 1e-4",
        "-snes_max_it 200",
        "-snes_PicardSwitchToNewton_rtol 5e-2",
        "-snes_NewtonSwitchToPicard_it 20",
        "-js_ksp_type fgmres",
        "-js_ksp_max_it 20",
        "-js_ksp_atol 1e-5",
        "-js_ksp_rtol 1e-4",
        "-snes_linesearch_type l2",
        "-snes_linesearch_maxstep 10",
        "-da_refine_y 1"
    ]
)

model.SolutionParams = SolutionParams(
    eta_min=1e-1 * refViscosity,
    eta_ref=refViscosity,
    eta_max=2e3 * refViscosity,
    min_cohes=1e6,
    init_guess=0  # initial guess flag
    # DII=1e-18
)
model.Time = Time(time_end=2000.0,
    dt=0.001,
    dt_min=0.00001,
    dt_max=0.1,
    nstep_max=400,
    nstep_out=10
)
run_lamem(model, 10)
