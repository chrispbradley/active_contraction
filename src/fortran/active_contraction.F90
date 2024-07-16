!> Main program
PROGRAM ActiveContractionExample

  USE OpenCMISS

  IMPLICIT NONE

  !Test program parameters

  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP

  INTEGER(OC_Intg), PARAMETER :: ContextUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: RegionUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(OC_Intg), PARAMETER :: MeshUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: DecomposerUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(OC_Intg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(OC_Intg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(OC_Intg), PARAMETER :: FieldGPUserNumber=5 ! temp/test
  INTEGER(OC_Intg), PARAMETER :: IndependentFieldUserNumber=6
  INTEGER(OC_Intg), PARAMETER :: EquationsSetFieldUserNumber=7

  INTEGER(OC_Intg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(OC_Intg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(OC_Intg), PARAMETER :: LinearMeshComponentNumber=2


  INTEGER(OC_Intg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(OC_Intg), PARAMETER :: ProblemUserNumber=1

  REAL(OC_RP), PARAMETER :: START_TIME = 0.0, END_TIME = 10.0, DT = 1  ! ms

  LOGICAL, PARAMETER :: TEST_GAUSS_POINT_FIELD = .FALSE.
  LOGICAL  :: directory_exists = .FALSE.

  !Program types

  !Program variables

  INTEGER(OC_Intg), PARAMETER, DIMENSION(1:27) :: ROTATE_ELEM = [ 1, 2, 3,10,11,12,19,20,21, &
                                                                  &  4, 5, 6,13,14,15,22,23,24, &  ! swap xi2 and xi3 directions for fiber angle
                                                                  &  7, 8, 9,16,17,18,25,26,27 ]  ! swap xi2 and xi3 directions for fiber angle

  INTEGER(OC_Intg) :: DecompositionIndex,EquationsSetIndex  
  INTEGER(OC_Intg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(OC_Intg) :: D, E, N, I, num_args, ix

  REAL(OC_RP) :: TMP
  REAL(OC_RP), DIMENSION(1:7) :: COSTA_PARAMS =  [ 0.2, 30.0, 12.0, 14.0, 14.0, 10.0, 18.0 ] ! a bff bfs bfn bss bsn bnn

  INTEGER(OC_Intg), dimension(:,:), allocatable :: Elements
  REAL(OC_RP)     , dimension(:,:), allocatable :: Nodes
  REAL(OC_RP)     , dimension(:,:), allocatable :: DirichletConditions
  REAL(OC_RP)     , dimension(:,:), allocatable :: Fibers
  REAL(OC_RP)     , dimension(:,:), allocatable :: ActivationTimes

  !OpenCMISS variables
  TYPE(OC_BasisType) :: QuadraticBasis, LinearBasis
  TYPE(OC_BoundaryConditionsType) :: BoundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: CoordinateSystem
  TYPE(OC_MeshType) :: Mesh
  TYPE(OC_DecompositionType) :: Decomposition
  TYPE(OC_DecomposerType) :: Decomposer
  TYPE(OC_EquationsType) :: Equations
  TYPE(OC_EquationsSetType) :: EquationsSet
  TYPE(OC_FieldType) :: GeometricField,EquationsSetField,FibreField,MaterialField,DependentField, GPfield, IndependentField
  TYPE(OC_FieldsType) :: Fields
  TYPE(OC_ProblemType) :: Problem
  TYPE(OC_RegionType) :: Region,WorldRegion
  TYPE(OC_SolverType) :: Solver,LinearSolver
  TYPE(OC_SolverEquationsType) :: SolverEquations
  TYPE(OC_MeshElementsType) :: QuadraticElements,LinearElements
  TYPE(OC_NodesType) :: CMNodes
  TYPE(OC_ControlLoopType) :: ControlLoop
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  character(len=256), dimension(:), allocatable :: args

  !Generic OpenCMISS variables
  INTEGER(OC_Intg) :: Err

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(contextUserNumber,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  num_args = command_argument_count()
  allocate(args(num_args))

  DO ix = 1, num_args
    CALL get_command_argument(ix,args(ix))
  END DO

  IF (num_args .gt. 0) THEN
    open(unit = 2, file = args(1))
  ELSE
!   open(unit = 2, file = "./input/hollowcylq-221.in")
    open(unit = 2, file = "./lvq-842.in")
  ENDIF
  IF (num_args .gt. 1) THEN
    open(unit = 3, file = args(2))
  ELSE
!   open(unit = 3, file = "./input/hollowcylq-221.in.gpactiv")
    open(unit = 3, file = "./lvq-842.in.gpactiv")
  ENDIF
  CALL read_mesh(2, Elements, Nodes, DirichletConditions, Fibers)
  CALL read_activation_times(3, ActivationTimes)
  close(2)
  close(3)

  NumberOfDomains=NumberOfComputationalNodes


  !Create a 3D rectangular cartesian coordinate system
  CALL OC_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL OC_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,Err)
  CALL OC_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL OC_Region_Initialise(Region,Err)
  CALL OC_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL OC_Region_LabelSet(Region,"Region",Err)
  CALL OC_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL OC_Region_CreateFinish(Region,Err)

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL OC_Basis_Initialise(LinearBasis,Err)
  CALL OC_Basis_CreateStart(LinearBasisUserNumber,context,LinearBasis,Err)
  CALL OC_Basis_CreateFinish(LinearBasis,Err)

  CALL OC_Basis_Initialise(QuadraticBasis,Err)
  CALL OC_Basis_CreateStart(QuadraticBasisUserNumber,context,QuadraticBasis,Err)
  CALL OC_Basis_InterpolationXiSet(QuadraticBasis,[OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL OC_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [OC_BASIS_MID_QUADRATURE_SCHEME,OC_BASIS_MID_QUADRATURE_SCHEME,OC_BASIS_MID_QUADRATURE_SCHEME],Err)
  CALL OC_Basis_CreateFinish(QuadraticBasis,Err)

  !Create a mesh with two components, Quadratic for geometry and fibers and linear lagrange
  !for hydrostatic pressure and material properties
  CALL OC_Mesh_Initialise(Mesh,Err)
  CALL OC_Mesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err) ! dim = 3
  CALL OC_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  CALL OC_Mesh_NumberOfElementsSet(Mesh,size(Elements,2),Err) ! num elts
  
  !define nodes for the mesh
  CALL OC_Nodes_Initialise(CMNodes,Err)
  CALL OC_Nodes_CreateStart(Region,size(Nodes,2),CMNodes,Err) ! num nodes
  CALL OC_Nodes_CreateFinish(CMNodes,Err)
  !Quadratic component : from file
  CALL OC_MeshElements_Initialise(QuadraticElements,Err)
  CALL OC_MeshElements_CreateStart(Mesh,QuadraticMeshComponentNumber,QuadraticBasis,QuadraticElements,Err)
  DO E=1,size(Elements,2)
    CALL OC_MeshElements_NodesSet(QuadraticElements,E, Elements(ROTATE_ELEM,E),Err)
  ENDDO
  CALL OC_MeshElements_CreateFinish(QuadraticElements,Err)
  !linear Lagrange component: numbers do not need to be continuous from 1? -> use quadratic nodeno on corners
  CALL OC_MeshElements_Initialise(LinearElements,Err)
  CALL OC_MeshElements_CreateStart(Mesh,LinearMeshComponentNumber,LinearBasis,LinearElements,Err)
  DO E=1,size(Elements,2)
    CALL OC_MeshElements_NodesSet(LinearElements,E, Elements( [1,3,7,9,19,21,25,27],E),Err)
  ENDDO
  CALL OC_MeshElements_CreateFinish(LinearElements,Err)

  !finish mesh creation
  CALL OC_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(Decomposition,Err)
  CALL OC_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL OC_Decomposition_CreateFinish(Decomposition,Err)

  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(decomposerUserNumber,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL OC_Field_Initialise(GeometricField,Err)
  CALL OC_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL OC_Field_DecompositionSet(GeometricField,Decomposition,Err)
  CALL OC_Field_TypeSet(GeometricField,OC_FIELD_GEOMETRIC_TYPE,Err)  
  CALL OC_Field_NumberOfVariablesSet(GeometricField,1,Err) ! 1 var
  CALL OC_Field_NumberOfComponentsSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,3,Err)   ! 3 components of geom field
  CALL OC_Field_VariableLabelSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL OC_Field_ComponentMeshComponentSet(GeometricField,OC_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL OC_Field_CreateFinish(GeometricField,Err)

  !Set node positions
  DO N=1,size(Nodes,2)
  DO D=1,3
    CALL OC_Field_ParameterSetUpdateNode(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,N,D,Nodes(D, &
      & N),Err)
  ENDDO
  ENDDO

  !Create a fibre field and attach it to the geometric field  
  CALL OC_Field_Initialise(FibreField,Err)
  CALL OC_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL OC_Field_TypeSet(FibreField,OC_FIELD_FIBRE_TYPE,Err)
  CALL OC_Field_DecompositionSet(FibreField,Decomposition,Err)        
  CALL OC_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL OC_Field_NumberOfComponentsSet(FibreField,OC_FIELD_U_VARIABLE_TYPE,3,Err)   ! 1 var, 3 components -> angles!
  CALL OC_Field_VariableLabelSet(FibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  DO D=1,3
    CALL OC_Field_ComponentMeshComponentSet(FibreField,OC_FIELD_U_VARIABLE_TYPE,D,QuadraticMeshComponentNumber,Err) ! quadratic interp
    CALL OC_Field_ComponentInterpolationSet(FibreField,OC_FIELD_U_VARIABLE_TYPE,D,OC_FIELD_NODE_BASED_INTERPOLATION,Err) ! node based
  ENDDO
  CALL OC_Field_CreateFinish(FibreField,Err)

  !Set fiber directions
  DO N=1,size(Nodes,2)
  DO D=1,3
    Fibers(D,N) = 0 ! input file gives unit vectors and opencmiss expects angles. TODO: fix.
    CALL OC_Field_ParameterSetUpdateNode(FibreField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,1,N,D,Fibers(D,N), &
      & Err)
  ENDDO
  ENDDO

 ! create the gauss point based field, for testing
  IF(TEST_GAUSS_POINT_FIELD) THEN
  WRITE(*,*) '---------<TESTING GAUSS POINT FIELD>---------'
  CALL OC_Field_Initialise(GPfield,Err)
  CALL OC_Field_CreateStart(FieldGPUserNumber,Region,GPfield,Err)
  CALL OC_Field_TypeSet(GPfield,OC_FIELD_GENERAL_TYPE,Err) ! ?
  CALL OC_Field_DecompositionSet(GPfield,Decomposition,Err)        
  CALL OC_Field_GeometricFieldSet(GPfield,GeometricField,Err)

  CALL OC_Field_NumberOfComponentsSet(GPfield,OC_FIELD_U_VARIABLE_TYPE,2,Err)

  ! FOR COMPONENTS
  CALL OC_Field_ComponentInterpolationSet(GPfield,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! GP based
  CALL OC_Field_ComponentInterpolationSet(GPfield,OC_FIELD_U_VARIABLE_TYPE,2,OC_FIELD_GAUSS_POINT_BASED_INTERPOLATION,Err) ! GP based

  CALL OC_Field_CreateFinish(GPfield,Err)

  CALL OC_Field_ComponentValuesInitialise(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,3.14_OC_RP,Err) ! init!
  CALL OC_Field_ComponentValuesInitialise(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,2.17_OC_RP,Err) ! init!
  CALL OC_Field_ComponentValuesInitialise(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,4.14_OC_RP,Err) ! set to const

  ! test gauss point field
  D=0;
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL OC_Field_ParameterSetGetGaussPoint(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I,E,1,TMP,Err)
    CALL OC_Field_ParameterSetUpdateGaussPoint(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I,E,1,TMP+D,Err)
    CALL OC_Field_ParameterSetGetGaussPoint(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I,E,2,TMP,Err)
    CALL OC_Field_ParameterSetUpdateGaussPoint(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I,E,2,TMP+D,Err)
    D=D+1
  ENDDO
  ENDDO

  D=0;
  DO E=1,size(Elements,2)
  DO I=1,8
    CALL OC_Field_ParameterSetGetGaussPoint(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I,E,1,TMP,Err)
    WRITE(*,*) 'COMPONENT 1 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 3.14 + 1 + D
    CALL OC_Field_ParameterSetGetGaussPoint(GPfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I,E,2,TMP,Err)
    WRITE(*,*) 'COMPONENT 2 ELEMENT ',E,', GP ', I, ' = ',TMP, ' Should = ', 2.17 + D
    D=D+1;
  ENDDO
  ENDDO
  WRITE(*,*) '---------</TESTING GAUSS POINT FIELD>---------'
  ENDIF  ! TEST GP FIELD

  !Create the equations_set
  CALL OC_Field_Initialise(EquationsSetField,Err)
  CALL OC_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  ! CHANGED
  CALL OC_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL OC_Field_Initialise(DependentField,Err)
  CALL OC_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL OC_Field_VariableLabelSet(DependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL OC_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL OC_Field_Initialise(MaterialField,Err)
  CALL OC_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL OC_Field_VariableLabelSet(MaterialField,OC_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL OC_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Costa material parameters
  DO I=1,7
   CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,I, &
     & COSTA_PARAMS(I),Err)
  END DO

!  CALL OC_Field_ComponentValuesInitialise(MaterialField,OC_FIELD_V_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2.0_OC_RP,Err) ! activate at time 2. TODO: inhomogeneous
  ! inhomogeneous activation times from file
  DO E=1,size(Elements,2)
  DO I=1,27
    CALL OC_Field_ParameterSetUpdateGaussPoint(MaterialField,OC_FIELD_V_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & I, E, 1, ActivationTimes(E,ROTATE_ELEM(I)), Err) ! rotating an element with 27 nodes can be done in the same way as the Gauss points
  ENDDO
  ENDDO

  ! create independent field
  CALL OC_Field_Initialise(IndependentField,Err)
  CALL OC_EquationsSet_IndependentCreateStart(EquationsSet,IndependentFieldUserNumber,IndependentField,Err)
  CALL OC_EquationsSet_IndependentCreateFinish(EquationsSet,Err)


  !Create the equations set equations
  CALL OC_Equations_Initialise(Equations,Err)
  CALL OC_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL OC_Equations_SparsityTypeSet(Equations,OC_EQUATIONS_SPARSE_MATRICES,Err)
  CALL OC_Equations_OutputTypeSet(Equations,OC_EQUATIONS_NO_OUTPUT,Err)
  CALL OC_EquationsSet_EquationsCreateFinish(EquationsSet,Err)   

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,Err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,Err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(GeometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,Err)
  CALL OC_Field_ComponentValuesInitialise(DependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4,-0.0_OC_RP, &
    & Err) ! -8?

  !Define the problem
  CALL OC_Problem_Initialise(Problem,Err)
  CALL OC_Problem_CreateStart(ProblemUserNumber,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE],Problem,Err)
   ! CHANGED TO OC_PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE
  CALL OC_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL OC_ControlLoop_Initialise(ControlLoop,Err)
  CALL OC_Problem_ControlLoopCreateStart(Problem,Err)
   CALL OC_Problem_ControlLoopGet(Problem,OC_CONTROL_LOOP_NODE,ControlLoop,Err)

   CALL OC_ControlLoop_TimesSet(ControlLoop, START_TIME - DT, END_TIME, DT, Err) ! set begin/end timings  . START AT -DT TO SOLVE FOR 0 AS WELL
   CALL OC_ControlLoop_TimeOutputSet(ControlLoop,1,Err)    !Set the output timing
  CALL OC_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_Solver_Initialise(LinearSolver,Err)
  CALL OC_Problem_SolversCreateStart(Problem,Err)
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL OC_Solver_OutputTypeSet(Solver,OC_SOLVER_PROGRESS_OUTPUT,Err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(Solver,OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err) ! faster than OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED ?
  CALL OC_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL OC_Solver_LinearTypeSet(LinearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL OC_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(Solver,Err)
  CALL OC_SolverEquations_Initialise(SolverEquations,Err)
  CALL OC_Problem_SolverEquationsCreateStart(Problem,Err)   
  CALL OC_Problem_SolverGet(Problem,OC_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL OC_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL OC_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL OC_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)


  DO I=1,size(DirichletConditions,2)
    N = INT(DirichletConditions(1,I))
    D = INT(DirichletConditions(2,I))
    CALL OC_BoundaryConditions_SetNode(BoundaryConditions,DependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,N,D,&
         & OC_BOUNDARY_CONDITION_FIXED, Nodes(D,N) + DirichletConditions(3,I),Err)  ! current + offset
  ENDDO


  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL OC_Problem_Solve(Problem,Err)

  INQUIRE(file="./results", exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir ./results")
  END IF

  !Output solution  
  CALL OC_Fields_Initialise(Fields,Err)
  CALL OC_Fields_Create(Region,Fields,Err)
  CALL OC_Fields_NodesExport(Fields,"./results/ActiveContraction","FORTRAN",Err)
  CALL OC_Fields_ElementsExport(Fields,"./results/ActiveContraction","FORTRAN",Err)
  CALL OC_Fields_Finalise(Fields,Err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP


contains
  subroutine read_mesh(fp, elements, node_coords, fixed_nodes, fibers)
    INTEGER(OC_Intg), intent(in)    :: fp  !< file 'pointer'
    INTEGER(OC_Intg), dimension(:,:), allocatable, intent(inout) :: elements !< element topology
    REAL(OC_RP)     , dimension(:,:), allocatable, intent(inout) :: node_coords  !< initial positions etc
    REAL(OC_RP)     , dimension(:,:), allocatable, intent(inout) :: fixed_nodes  !< dirichlet boundary conditions
    REAL(OC_RP)     , dimension(:,:), allocatable, intent(inout) :: fibers       !< unit vectors for fiber dir

    INTEGER(OC_Intg) :: number_of_elts, number_of_dims, number_of_fixednodes, maxnodenr, i,j, eltno, nodeno

    character(len=256) :: desc_str
 
    maxnodenr = 0

    read (fp,*) desc_str
    write (*,*) "'number_of_elements':", desc_str

    read (fp,*) number_of_elts
    write (*,*) "# of elements:", number_of_elts

    allocate(elements(1:27,1:number_of_elts))

    read (fp,*) desc_str
    write (*,*) "'topology':", desc_str
 
    do i = 1, number_of_elts
      read (fp,*)  eltno, elements(:,i)
      write (*,*) "nodes of element ", eltno, " are : ", elements(:,i)
      do j=1,27
        maxnodenr = max(maxnodenr,elements(j,i))
      end do
    end do
   

    read (fp,*) desc_str
    read (fp,*) number_of_dims
    write (*,*) "max node nr:", maxnodenr, " dimensions: ", number_of_dims
    write (*,*) "'initial_positions':", desc_str
 
    allocate(node_coords(1:number_of_dims,1:maxnodenr))
    do i = 1, maxnodenr
      read (fp,*) nodeno, node_coords(:,i)  ! expect input in order
      write (*,*) "initial value of node ", nodeno, "=", i ," is ", node_coords(:,i) 
    end do

    read (fp,*) desc_str
    write (*,*) "'fixed_nodes':", desc_str
    read (fp,*) number_of_fixednodes
    write (*,*) "# of dirichlet bc:", number_of_fixednodes
    allocate(fixed_nodes(3,number_of_fixednodes))
    do i = 1, number_of_fixednodes
      read (fp,*) fixed_nodes(:,i)
      write (*,*) "dirichlet bc: ", fixed_nodes(:,i) 
    end do
  
    read (fp,*) desc_str ! traction
    read (fp,*) desc_str ! 0 0
    read (fp,*) desc_str
    write (*,*) "'fibers':", desc_str

    allocate(fibers(3,1:maxnodenr))
    do i = 1, maxnodenr
      read (fp,*) nodeno, fibers(:,i)  ! expect input in order
      write (*,*) "fiber dir at ", nodeno, "=", i ," is ", fibers(:,i) 
    end do


  end subroutine read_mesh    


  subroutine read_activation_times(fp, activtime)
    INTEGER(OC_Intg), intent(in)    :: fp  !< file 'pointer'
    REAL(OC_RP), dimension(:,:), allocatable, intent(inout) :: activtime !< elements x 27 array of activation times
    INTEGER(OC_Intg)   :: num_elt, num_gp, I
    read (fp,*) num_elt, num_gp
    write(*,*) 'reading activation times on ',num_elt,' elements / ',num_gp,' Gauss points'
    allocate( activtime(1:num_elt,1:num_gp))
    do I=1,num_elt
      read (fp,*) activtime(I,:)
      write(*,*) 'element ',I,' : ',activtime(I,:)
    enddo
  end subroutine read_activation_times

END PROGRAM ActiveContractionExample

