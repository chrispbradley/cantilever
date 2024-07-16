#> Main script
import os

# Intialise OpenCMISS-Iron
from opencmiss.opencmiss import OpenCMISS_Python as oc

# Set problem parameters

width = 60.0
length = 40.0
height = 40.0
density=9.0E-4 #in g mm^-3
gravity=[0.0,0.0,-9.81] #in m s^-2

usePressureBasis = False
NumberOfGaussXi = 2

contextUserNumber = 1
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
decomposerUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
materialFieldUserNumber = 3
dependentFieldUserNumber = 4
sourceFieldUserNumber = 5
equationsSetFieldUserNumber = 6
equationsSetUserNumber = 1
problemUserNumber = 1

context = oc.Context()
context.Create(contextUserNumber)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

# Set all diganostic levels on for testing
#oc.DiagnosticsSetOn(oc.DiagnosticTypes.All,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

numberOfLoadIncrements = 1
numberGlobalXElements = 1
numberGlobalYElements = 1
numberGlobalZElements = 1
InterpolationType = 1
if(numberGlobalZElements==0):
    numberOfXi = 2
else:
    numberOfXi = 3

# Get the number of computational nodes and this computational node number
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,context)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = oc.Region()
region.CreateStart(regionUserNumber,worldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define basis
basis = oc.Basis()
basis.CreateStart(basisUserNumber,context)
if InterpolationType in (1,2,3,4):
    basis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
elif InterpolationType in (7,8,9):
    basis.type = oc.BasisTypes.SIMPLEX
basis.numberOfXi = numberOfXi
basis.interpolationXi = [oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
if(NumberOfGaussXi>0):
    basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
basis.CreateFinish()

if(usePressureBasis):
    # Define pressure basis
    pressureBasis = oc.Basis()
    pressureBasis.CreateStart(pressureBasisUserNumber,context)
    if InterpolationType in (1,2,3,4):
        pressureBasis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
    elif InterpolationType in (7,8,9):
        pressureBasis.type = oc.BasisTypes.SIMPLEX
    pressureBasis.numberOfXi = numberOfXi
    pressureBasis.interpolationXi = [oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
    if(NumberOfGaussXi>0):
        pressureBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
    pressureBasis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = oc.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = oc.GeneratedMeshTypes.REGULAR
if(usePressureBasis):
    generatedMesh.basis = [basis,pressureBasis]
else:
    generatedMesh.basis = [basis]
if(numberGlobalZElements==0):
    generatedMesh.extent = [width,height]
    generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements]
else:
    generatedMesh.extent = [width,length,height]
    generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]
# Finish the creation of a generated mesh in the region
mesh = oc.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = oc.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.CreateFinish()

# Decompose 
decomposer = oc.Decomposer()
decomposer.CreateStart(decomposerUserNumber,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

# Create a field for the geometry
geometricField = oc.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(oc.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(oc.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
if InterpolationType == 4:
    geometricField.fieldScalingType = oc.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = oc.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(oc.FieldTypes.FIBRE)
fibreField.DecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(oc.FieldVariableTypes.U,"Fibre")
if InterpolationType == 4:
    fibreField.fieldScalingType = oc.FieldScalingTypes.ARITHMETIC_MEAN
fibreField.CreateFinish()

# Create the equations_set
equationsSetField = oc.Field()
equationsSet = oc.EquationsSet()
equationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY,
    oc.EquationsSetTypes.FINITE_ELASTICITY,
    oc.EquationsSetSubtypes.MOONEY_RIVLIN]
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create the dependent field
dependentField = oc.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(oc.FieldVariableTypes.U,"Dependent")
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U,4,oc.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.T,4,oc.FieldInterpolationTypes.ELEMENT_BASED)
if(usePressureBasis):
    # Set the pressure to be nodally based and use the second mesh component
    if InterpolationType == 4:
        dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U,4,oc.FieldInterpolationTypes.NODE_BASED)
        dependentField.ComponentInterpolationSet(oc.FieldVariableTypes.T,4,oc.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,4,2)
    dependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.T,4,2)
if InterpolationType == 4:
    dependentField.fieldScalingType = oc.FieldScalingTypes.ARITHMETIC_MEAN
equationsSet.DependentCreateFinish()


# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1)
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2)
oc.Field.ParametersToFieldParametersComponentCopy(
    geometricField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3,
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3)
oc.Field.ComponentValuesInitialiseDP(
    dependentField,oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,4,0.0)

# Create the material field
materialField = oc.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(oc.FieldVariableTypes.U,"Material")
materialField.VariableLabelSet(oc.FieldVariableTypes.V,"Density")
equationsSet.MaterialsCreateFinish()

# Set Mooney-Rivlin constants c10 and c01 respectively.
materialField.ComponentValuesInitialiseDP(
    oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,2.0)
materialField.ComponentValuesInitialiseDP(
    oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,2.0)
materialField.ComponentValuesInitialiseDP(
    oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,1,density)

#Create the source field with the gravity vector
sourceField = oc.Field()
equationsSet.SourceCreateStart(sourceFieldUserNumber,sourceField)
if InterpolationType == 4:
    sourceField.fieldScalingType = oc.FieldScalingTypes.ARITHMETIC_MEAN
else:
    sourceField.fieldScalingType = oc.FieldScalingTypes.UNIT
equationsSet.SourceCreateFinish()

#Set the gravity vector component values
sourceField.ComponentValuesInitialiseDP(
    oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,gravity[0])
sourceField.ComponentValuesInitialiseDP(
    oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,gravity[1])
sourceField.ComponentValuesInitialiseDP(
    oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3,gravity[2])

# Create equations
equations = oc.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = oc.EquationsSparsityTypes.SPARSE
equations.outputType = oc.EquationsOutputTypes.NONE
equations.outputType = oc.EquationsOutputTypes.ELEMENT_MATRIX
equationsSet.EquationsCreateFinish()

# Define the problem
problem = oc.Problem()
problemSpecification = [oc.ProblemClasses.ELASTICITY,
                        oc.ProblemTypes.FINITE_ELASTICITY,
                        oc.ProblemSubtypes.STATIC_FINITE_ELASTICITY]
problem.CreateStart(problemUserNumber,context,problemSpecification)
problem.CreateFinish()

# Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = oc.ControlLoop()
problem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
problem.ControlLoopCreateFinish()

# Create problem solver
nonLinearSolver = oc.Solver()
linearSolver = oc.Solver()
problem.SolversCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,nonLinearSolver)
nonLinearSolver.outputType = oc.SolverOutputTypes.MONITOR
nonLinearSolver.outputType = oc.SolverOutputTypes.MATRIX
nonLinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
nonLinearSolver.NewtonSolutionToleranceSet(1e-14)
nonLinearSolver.NewtonRelativeToleranceSet(1e-14)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = oc.LinearSolverTypes.DIRECT
#linearSolver.libraryType = oc.SolverLibraries.LAPACK
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = oc.Solver()
solverEquations = oc.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = oc.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = oc.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
# Set x=0 nodes to no x displacment
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,1,1,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,3,1,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,5,1,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,7,1,oc.BoundaryConditionsTypes.FIXED,0.0)

# Set y=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,1,2,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,3,2,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,5,2,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,7,2,oc.BoundaryConditionsTypes.FIXED,0.0)

# Set z=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,1,3,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,3,3,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,5,3,oc.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,7,3,oc.BoundaryConditionsTypes.FIXED,0.0)

# Set z=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,6,3,oc.BoundaryConditionsTypes.FIXED,-5.0)
boundaryConditions.AddNode(dependentField,oc.FieldVariableTypes.U,1,1,8,3,oc.BoundaryConditionsTypes.FIXED,-5.0)

solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

if not os.path.exists("./results"):
    os.makedirs("./results")

# Export results
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/Cantilever","FORTRAN")
fields.ElementsExport("./results/Cantilever","FORTRAN")
fields.Finalise()

