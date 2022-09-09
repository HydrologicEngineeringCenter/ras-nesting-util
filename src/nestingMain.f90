!*********************
program nestingMain
!*********************
  use nestingModule
  implicit none
  character(len=256):: parentGridFile
  character(len=256):: childGridFile
  character(len=256):: pointsCoordFile !Points coordinate csv file and path
  character(len=256):: polyfile        !Bounding polygon csv file and path
  character(len=256):: savePointFile
  character(len=256),allocatable:: parentWSEFiles(:)
  character(len=256),allocatable:: childWSEFiles(:)
  character(len=256),allocatable:: pointsWSEFiles(:)
  character(len=256),allocatable:: arg(:)
  character(len=256):: msg
  integer:: status, narg, i, nts
  double precision:: factor
  
  narg = command_argument_count()
  if(narg == 0)then
    call TEST_triangulate()
    call TEST_INPolygon()
    !call TEST_gensub()
    !call TEST_extractwse()
    !call TEST_interpwse()
    call TEST_trisp()
    !call TEST_extractsp()
    !call TEST_interpsp()
    stop
  endif
  
  allocate(arg(narg))
  
  do i=1,narg
    call get_command_argument(i,arg(i))
  enddo
  
  !Command 
  selectcase(arg(1))
  case('-gensub')
    if(narg /= 4)then
      write(*,*) 'ERROR: Invalid number of arguments'
      pause
      stop
    endif
    parentGridFile = arg(2)
    polyFile = arg(3)
    childGridFile = arg(4)
    call nesting_gensub(parentGridFile, &
      polyFile, childGridFile, status, msg)
    
  case('-extractwse')
    if(narg < 4)then
      write(*,*) 'ERROR: Invalid number of arguments'
      pause
      stop
    endif
    nts = (narg - 2)/2
    childGridFile = arg(2)
    allocate(parentWSEFiles(nts),childWSEFiles(nts))
    do i=1,nts
      parentWSEFiles(i) = arg(2+i)
      childWSEFiles(i) = arg(3+i)
    enddo
    call nesting_extractwse(childGridFile, &
      nts, parentWSEFiles, childWSEFiles, status, msg)
    
  case('-interpwse')
    if(narg < 5)then
      write(*,*) 'ERROR: Invalid number of arguments'
      pause
      stop
    endif
    nts = (narg - 4)/2
    childGridFile = arg(2)
    pointsCoordFile = arg(3)
    allocate(childWSEFiles(nts),pointsWSEFiles(nts))
    do i=1,nts
      childWSEFiles(1) = arg(3+i)
      pointsWSEFiles(1) = arg(4+i)
    enddo
    read(arg(narg),*) factor
    call nesting_interpwse(childGridFile, pointsCoordFile, &
      nts, childWSEFiles, pointsWSEFiles, factor, status, msg)
    
  case('-trisp')
    if(narg /= 4)then
      write(*,*) 'ERROR: Invalid number of arguments'
      pause
      stop
    endif
    savePointFile = arg(2)
    polyFile = arg(3)
    childGridFile = arg(4)
    call nesting_trisp(savePointFile,polyFile,childGridFile,status,msg)
    
  endselect
  
  !Clean
  if(allocated(parentWSEFiles))then
    deallocate(parentWSEFiles)
  endif
  if(allocated(childWSEFiles))then
    deallocate(childWSEFiles)
  endif
  if(allocated(pointsWSEFiles))then
    deallocate(pointsWSEFiles)
  endif
  if(allocated(arg))then
    deallocate(arg)
  endif
  
  !Finish
  write(*,*)'Complete'
  
  stop
  
contains
  
!*****************************************************
  subroutine TEST_INPolygon()
!*****************************************************
    use nestingModule
    implicit none
    integer:: np
    double precision,allocatable:: xp(:),yp(:)
    double precision:: xi,yi
    integer:: status
    logical:: inPoly
    
    status = 0
    
    !Test 1
    np = 5
    allocate(xp(np),yp(np))
    xp = (/3.365826, 3.401575, 3.464305, 3.450711, 3.41/)
    yp = (/0.1425249, 0.1108501, 0.1149755, 0.1455467, 0.17/)
    xi = 3.400732
    yi = 0.1459093
    inPoly = inPolygon(np,xp,yp,xi,yi)
    if(.not.inPoly)then
      status = -1
    endif
    
    !Test 3
    xp = (/3.365826, 3.401575, 3.464305, 3.450711/)
    yp = (/0.1425249, 0.1108501, 0.1149755, 0.1455467/)
    xi = 3.398772
    yi = 0.1436978
    inPoly = inPolygon(np,xp,yp,xi,yi)
    if(inPoly)then
      status = -2
    endif
    
    !Clean
    deallocate(xp,yp)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_INPolygon'
    endif
    
  endsubroutine
  
!****************************************
  subroutine TEST_gensub()
!****************************************
    use nestingModule
    implicit none
    integer:: status
    character(len=256):: msg
    character(len=512):: parentGridFile,polyFile,childGridFile
  
    parentGridFile = "C:\Projects\hec-ras\engine\Nesting\RunTime\fort.14"
    polyFile = "C:\Projects\hec-ras\engine\Nesting\RunTime\box.csv"
    childGridFile = "C:\Projects\hec-ras\engine\Nesting\RunTime\ChildBox.grd"
    
    !parentGridFile = "C:\Projects\nesting\RunTime\S2G_new_lidar_full_roads_PAOF_NoProjects_v7.grd"
    !polyFile = "C:\Projects\nesting\RunTime\polygonNeches.csv"
    !childGridFile = "C:\Projects\nesting\RunTime\ChildNeches.grd"
    
    !parentGridFile = "D:\HEC\HEC-RAS\Shahidul\RAS_enhancement\Mesh_Allflows\CCDB_all_CC_N_DB_v2.grd"
    !polyFile = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\polygonBoundaryNeches.csv"
    !childGridFile = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\ChildNeches.grd"
    
    call nesting_gensub(parentGridFile,polyFile,childGridFile,status,msg)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_gensub: ', trim(msg)
    endif
    
  endsubroutine
  
!****************************************
  subroutine TEST_trisp()
! Triangulate Save Points
!****************************************
    use nestingModule
    implicit none
    integer:: status
    character(len=256):: msg
    character(len=512):: savePointFile,polyFile,childGridFile
    
    !savePointFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Savepoints\elev_stat.151"
    !polyFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\boxBoundaryNeches.csv"
    !childGridFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\SubdomainSavePointsNeches.grd"
    
    savePointFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Savepoints\elev_stat.151"
    polyFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\DickinsonBayou\Ike\boxBC.csv"
    childGridFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\DickinsonBayou\Ike\DickinsonChildIkeSavePoints.grd"
    
    call nesting_trisp(savePointFile,polyFile,childGridFile,status,msg)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_gensub: ', trim(msg)
    endif
    
  endsubroutine
  
!****************************************
  subroutine TEST_extractsp()
! Extract water surface elevations from 
! parent or child fort.63 file
!****************************************
    use nestingModule
    implicit none
    integer:: status
    character(len=256):: msg
    character(len=512):: childGridFile,parentWSEFiles(1),childWSEFiles(1)
  
    !childGridFile = "C:\Projects\nesting\RunTime\Child.grd"
    !parentWSEFiles(1) = "C:\Projects\nesting\RunTime\fort.63"
    !childWSEFiles(1) = "C:\Projects\nesting\RunTime\Child.63"
    
    !childGridFile = "C:\Projects\nesting\RunTime\ChildNeches.grd"
    !parentWSEFiles(1) = "C:\Projects\nesting\RunTime\CTXCS_TP_0011_HIS_Tides_1_SLC_0_RFC_0_WAV_1_GCP_S2G03BE01_fort.63"
    !childWSEFiles(1) = "C:\Projects\nesting\RunTime\ChildNeches.63"
    
    childGridFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\SubdomainSavePointsNeches.grd"
    parentWSEFiles(1) = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Savepoints\" &
      // "CTXCS_TP_0009_HIS_Tides_1_SLC_0_RFC_0_WAV_1_GCP_CCDB01E02T_All_Flows_CCDB_TIMEVARY_Oct13_fort.61"
    childWSEFiles(1) = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\Ike\ChildNechesIke.61"
    
    call nesting_extractwse(childGridFile,1,parentWSEFiles,childWSEFiles,status,msg)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_extractwse: ', trim(msg)
    endif
    
  endsubroutine

!****************************************
  subroutine TEST_extractwse()
! Extract water surface elevations from 
! parent or child fort.63 file
!****************************************
    use nestingModule
    implicit none
    integer:: status
    character(len=256):: msg
    character(len=512):: childGridFile,parentWSEFiles(1),childWSEFiles(1)
  
    !childGridFile = "C:\Projects\nesting\RunTime\Child.grd"
    !parentWSEFiles(1) = "C:\Projects\nesting\RunTime\fort.63"
    !childWSEFiles(1) = "C:\Projects\nesting\RunTime\Child.63"
    
    !childGridFile = "C:\Projects\nesting\RunTime\ChildNeches.grd"
    !parentWSEFiles(1) = "C:\Projects\nesting\RunTime\CTXCS_TP_0011_HIS_Tides_1_SLC_0_RFC_0_WAV_1_GCP_S2G03BE01_fort.63"
    !childWSEFiles(1) = "C:\Projects\nesting\RunTime\ChildNeches.63"
    
    childGridFile = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\ChildNeches.grd"
    parentWSEFiles(1) = "D:\HEC\HEC-RAS\Shahidul\RAS_enhancement\Model_results\" &
      // "CTXCS_TP_0009_HIS_Tides_1_SLC_0_RFC_0_WAV_1_GCP_CCDB01E02T_All_Flows_CCDB_TIMEVARY_Oct13_fort.63"
    childWSEFiles(1) = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\Ike\ChildNechesIke.63"
    
    call nesting_extractwse(childGridFile,1,parentWSEFiles,childWSEFiles,status,msg)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_extractwse: ', trim(msg)
    endif
    
  endsubroutine
  
!****************************************
  subroutine TEST_interpwse()
!****************************************
    use nestingModule
    implicit none
    integer:: status
    character(len=256):: msg
    character(len=512):: childGridFile,pointsCoordFile
    character(len=512):: childWSEFiles(1),pointsWSEFiles(1)
    double precision:: factor
    
    !childGridFile = "C:\Projects\nesting\RunTime\Child.grd"
    !pointsCoordFile = "C:\Projects\nesting\RunTime\Boundary.csv"
    !childWSEFiles(1) = "C:\Projects\nesting\RunTime\Child.63"
    !pointsWSEFiles(1) = "C:\Projects\nesting\RunTime\Child.csv"
    
    !childGridFile = "C:\Projects\nesting\RunTime\ChildNeches.grd"
    !pointsCoordFile = "C:\Projects\nesting\RunTime\boundaryCoordinatesNeches.csv"
    !childWSEFiles(1) = "C:\Projects\nesting\RunTime\ChildNeches.63"
    !pointsWSEFiles(1) = "C:\Projects\nesting\RunTime\wseBoundaryNeches.csv"
    
    childGridFile = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\ChildNeches.grd"
    pointsCoordFile = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\boundaryCoordinatesNeches.csv"
    childWSEFiles(1) = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\Ike\ChildNechesIke.63"
    pointsWSEFiles(1) = "D:\HEC\HEC-RAS\Shahidul\Neches\Nesting\Ike\wseBoundaryNechesIke.csv"
    factor = 1d0/0.3048d0 !Convert from m to ft
    
    call nesting_interpwse(childGridFile,pointsCoordFile,&
      1,childWSEFiles,pointsWSEFiles,factor,status,msg)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_extractwse: ', trim(msg)
    endif
    
  endsubroutine
  
!****************************************
  subroutine TEST_interpsp()
!****************************************
    use nestingModule
    implicit none
    integer:: status
    character(len=256):: msg
    character(len=512):: childGridFile,pointsCoordFile
    character(len=512):: childWSEFiles(1),pointsWSEFiles(1)
    double precision:: factor
    
    childGridFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\SubdomainSavePointsNeches.grd"
    pointsCoordFile = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\boundaryCoordinatesNeches.csv"
    childWSEFiles(1) = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\Ike\ChildNechesIke.61"
    pointsWSEFiles(1) = "D:\HEC\HEC-RAS\Projects\Variable_WSE_BC\Nesting\Ike\wseBoundaryNechesIkeSP.csv"
    factor = 1d0/0.3048d0 !Convert from m to ft
    
    call nesting_interpwse(childGridFile,pointsCoordFile,&
      1,childWSEFiles,pointsWSEFiles,factor,status,msg)
    
    if(status /= 0)then
      write(*,*) 'ERROR: Failed TEST_extractwse: ', trim(msg)
    endif
    
  endsubroutine
  
!******************************************
  subroutine TEST_triangulate()
!******************************************
    use delaunay2d, only: dtris2
    implicit none
    integer,parameter:: point_num = 5
    double precision:: point_xy(2,point_num)
    integer:: tri_num, num
    integer:: tri_vert(3,point_num*2),tri_vert_true(3,4)
    integer:: tri_nabe(3,point_num*2)
    
    point_xy(1,1:point_num) = (/1d0, 3d0, 5d0, 6d0, 4d0/)
    point_xy(2,1:point_num) = (/3d0, 2d0, 4d0, 3d0, 3d0/)
    
    call dtris2(point_num, point_xy, tri_num, tri_vert, tri_nabe)
    
    !Correct Connectivity List
    tri_vert_true(1:3,1) = (/1, 2, 5/)
    tri_vert_true(1:3,2) = (/1, 5, 3/)
    tri_vert_true(1:3,3) = (/5, 2, 4/)
    tri_vert_true(1:3,4) = (/3, 5, 4/)
    
    !Check for errors
    num = sum(tri_vert_true(1:3,1:4) - tri_vert(1:3,1:4))
    if(num > 0)then
      write(*,*) 'ERROR: Failed TEST_triangulate'
    endif
    
  endsubroutine
  
endprogram
  