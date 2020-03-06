program NVT_from_bulk
! Single component multiset fitting algorithm.
! Use: multiset-inverse-smarter.exe < multiset.in 
!----------------------------------------
! INPUT FILE FORMAT (muliset.in)
!----------------------------------------
! The 1st line of the input file should contain a single value: 
!   - the total number of data sets included in this data file that will 
!     subsequently be fitted to.
! The 2nd line of the input file should contain three values separarated 
! by spaces:
!   - The number iterations of fitting used to find the best fit equilibrium 
!     constants.
!   - The exponent used to adjust the equilibrium constant at each iteration. 
!     This should be a value less than 1.0 but greater than 0. It's recommended
!     that you start with an value of 0.5. If the fitting proceedure does not 
!     converge or the convergence criterion fluctuates wildly, try lowering the 
!     exponent value.
!   - monomer concentration used for the generation of the gibbs free energy of 
!     association.
! The 3rd line and onwards should contain the data sets being fitted to.
! Each data set should contain one line with three values separated by spaces:
!   - The total number of monomers in the simulation for that data set.
!   - The number of different cluster sizes that show up in the data set.
!   - The effective free volume in the simulation used to generate the data set.
! The subseqent lines for the data set should each contain three values 
! separated by a space:
!   - The cluster size.
!   - The average cluster frequency for that cluster size.
!   - The standard deviation in the average cluster frequency.
! 
! If you use this code, please cite the following paper:
!     L. A. Patel and J. T. Kindt. Cluster Free Energies from Simple 
!     Simulations of Small Numbers of Aggregants: Nucleation of Liquid 
!     MTBE from Vapor and Aqueous Phases. J. Chem. Theory Comput., 2017, 
!     13 (3), pp 1023-1033.
! 
!----------------------------------------
! DEFINING VARIABLE AND ARRAY TYPES
!----------------------------------------
implicit none ! Cancels default naming conventions
integer, parameter:: Nmax=160
integer, parameter:: Nsetmax=50
       ! to go beyond Nmax=170, will need to 
       ! do something about overflow of factorials

character(8) :: date
character(10) :: time
character(5) :: zone, iterString
character(100) :: startTime,finishTime
real*8 :: startTime_cpu, finishTime_cpu
! Write formats, e.g. write(*,FMT1) instead of write(*,*)
character(*), parameter:: FMT2 = "(I3,4(1x,E12.5E3),1x,I3,(1x,E12.5E2))"
character(*), parameter:: FMT3 = "(2(1x,I3),1(1x,E12.5E3),1(1x,F12.6))"
character(*), parameter:: FMT4 = "(2(1x,I3),1(1x,E12.5E3),2(1x,F12.6))"
character(*), parameter:: FMT5 = "(I3,8(1x,E12.5E3))"

real*8:: conc(Nmax,Nsetmax),currentconc(Nmax)

real*8:: aveN(Nmax,Nsetmax),aveNout(0:Nmax,Nsetmax)
real*8:: aveNconstrained(0:Nmax,Nsetmax)
real*8:: stdvN(Nmax,Nsetmax),NrelW(Nmax,Nsetmax),weightSum(Nmax)

real*8:: totconc(Nsetmax),V(Nsetmax) 
real*8:: Kapp(Nmax,Nsetmax),Kactual(Nmax)
real*8:: alpha,nc1 

integer:: iset
integer:: Nsets,nsamp(0:Nmax)
integer:: ni_in(Nsetmax),Nsim(Nsetmax),NsimMax
integer:: i,j,k,l,niter,nci_gt_0(Nsetmax)
integer:: ilist(0:Nmax+1,Nsetmax),revlist(0:Nmax+1,Nsetmax)
integer:: most_prob_comp(0:Nmax+1,Nsetmax)

! Keeps track of the monomers in a set that have a conc greater than 0.
logical:: ci_gt_0(Nmax,Nsetmax) 

100 format (20d8.1)

!----------------------------------------
! GETTING THE DATE AND TIME
call date_and_time(DATE=date,TIME=time,ZONE=zone)
startTime="# Date: "//date(1:4)//"."//date(5:6)//"."//date(7:8)//&
"\n # Time: "//time(1:2)//":"//time(3:4)//":"//time(5:10)//&
" (UTC"//zone(1:3)//":"//zone(4:5)//")"
call cpu_time(startTime_cpu)

nci_gt_0=0
ci_gt_0=.FALSE.
ilist=0
totconc=0.0d0
conc=0.0d0
Kapp=0.0d0
Kactual=0.0d0
aveNout=0.0d0
weightSum=0.0d0

!----------------------------------------
! READING DATA

write (*,*) "##########"
write (*,*) "# READING IN DATA"
write (*,*) "##########"

read (*,*) Nsets 
if (Nsets.gt. Nsetmax) then
  write (0,*) "WARNING: Total number of data sets (",Nsets,&
  ") exceeds the maximum (",Nsetmax,")."
  stop
endif
read (*,*) niter, alpha, nc1
! niter: max number of iterations to find best fit equilibrium constants
! alpha: exponent to adjust equilibrium constant at each iteration
! C1_out: monomer concentration for DeltaG print in output
write(*,*) "Monomer concentration for output is ",nc1," molec/nm^3."

do iset=1,nsets
  read (*,*) Nsim(iset),ni_in(iset),V(iset)
  write (*,*) Nsim(iset),ni_in(iset)
  ! Nsim: total number of monomers in the simulation
  ! ni_in: the total number of i-mers whose concentration from simulation 
  !    will be given
  ! V: volume
  if (NsimMax.lt.Nsim(iset)) then 
    NsimMax=Nsim(iset) 
  endif
  if (Nsim(iset).gt.Nmax) then
    write (0,*) "WARNING: Nsim (",Nsim,") is greater than Nmax (",Nmax,&
    ").\n    Adjust the Nmax parameter to equal Nsim."
    stop
  endif

  ! Reading in the monomer concentration first and asserting that it is a
  ! non-zero value.
  read (*,*) i,aveN(i,iset),stdvN(i,iset)
  NrelW(i,iset)=aveN(i,iset)/stdvN(i,iset)
  weightSum(i)=weightSum(i)+NrelW(i,iset)
  conc(i,iset)=aveN(i,iset)/V(iset)
  totconc(iset)=totconc(iset)+conc(1,iset)
  if (conc(1,iset).le.0.0d0) then
    write (0,*) "WARNING: Concentration of free monomer must be > 0.0d0.",&
    "\n   Current free monomer concentration = ",conc(1,iset)
    stop
  endif
  ci_gt_0(1,iset)=.TRUE. ! bool of whether conc > 0 for that i-mer
  nci_gt_0(iset)=1 ! number of i-mers with conc > 0

  ! ilist contains the list of all i-mers present; 
  ! ilist(nci_gt_0) is the biggest aggregate observed.
  ilist(0,iset)=0
  ilist(1,iset)=1   ! first in list is 1 for monomer
  Kapp(1,iset)=1.0d0

  ! Cycling through the rest of the i-mers
  do j=2,ni_in(iset)
    read (*,*) i,aveN(i,iset),stdvN(i,iset)
    NrelW(i,iset)=aveN(i,iset)/stdvN(i,iset)
    weightSum(i)=weightSum(i)+NrelW(i,iset)
    ! input bulk system concentration of i-mers
    conc(i,iset)=aveN(i,iset)/V(iset)
    if (conc(i,iset).gt.0.0d0) then 
      ! keep track of what aggregate sizes are present
      if (i.le.Nsim(iset)) then !! REDUNDANT?
      ! shouldnt this be true all the time by virtue of 
      ! the simulation cluster sizes, i, being unable to 
      ! exceed the total number of molecules in the 
      ! simulation, Nsim(iset)?
        ci_gt_0(i,iset)=.TRUE.
        nci_gt_0(iset)=nci_gt_0(iset)+1
        ilist(nci_gt_0(iset),iset)=i
      endif
    endif
    totconc(iset)=totconc(iset)+i*conc(i,iset)
    ! apparent equilibrium constant from simulation
    Kapp(i,iset)=conc(i,iset)/conc(1,iset)**i 
  enddo 
  ilist(nci_gt_0(iset)+1,iset)=Nsim(iset)+1
  ! make revlist: keeps track of index of equal or 
  ! next-lowest aggregate size present; for any 
  ! aggregate size i, equal or next-lowest aggregate size is
  ! ilist(revlist(i))
  revlist(:,iset)=nci_gt_0(iset)
  do j=0,Nmax+1
    do i=nci_gt_0(iset), 0, -1
      if (ilist(i,iset).gt.j) then
        revlist(j,iset)=i-1
      endif
    enddo
  enddo
  !write (*,*) "# total monomer concentration is:",totconc(iset)
  !write (*,*) "# volume is ",Nsim(iset),"/",totconc(iset)," = ",V(iset)
enddo  ! finish reading data sets

open(UNIT=11,FILE="MS-multiset-xKapp.xvg",form='FORMATTED',status='REPLACE',action='WRITE')
write (11,*) startTime
write (11,*) "## XVG file containing the Kapparent values as the law of mass"
write (11,*) "## action would predict based on the concentrations of the i-mers"
write (11,*) "## in each simulation."
write (11,*) "## COL | VARIABLE"
write (11,*) "##  1  | iset, the data set number."
write (11,*) "##  2  | i, the cluster size."
write (11,*) "##  3  | Kapparent, the equilibrium rate constant for the i-mer."
write (11,*) "##  4  | Delta_G, the change in free energy for the i-mer at a "
write (11,*) "##       monomer concentration of ",nc1," molecules/nm^3."
write (11,*) "##  5  | NrelW(i,iset)/weightSum(i), relative weighting of the "
write (11,*) "##       i-mer frequency:"
write (11,*) "##       ((stdvN(i,iset)/aveN(i,iset))**(-1))/weightSum(i)."
write (11,*) "##  "
do iset=1,nsets
  do i=1,Nsim(iset)
  write (11,FMT4) iset, i, Kapp(i,iset), -log(Kapp(i,iset))-(i-1)*log(nc1), &
  NrelW(i,iset)/weightSum(i)
  enddo
enddo
close(unit=11)

write (*,*) "##########"
write (*,*) "# FINISHED READING IN DATA"
write (*,*) "##########"

!----------------------------------------
! GENERATE Kactual USING A WEIGHTED GEOMETRIC MEAN OF Kapp
! This is done using a weighted geometric mean of the Kapp values that have
! been read in.
! <x>g = (Sum_i^n x_i^(w_i))^(1/(Sum_i^n w_i))
! where x_i are the i data points and w_i are the respective weightings which in
! this case we have based off of the standard deviation of the individual Kapp
! values which in turn is calculated from the standard devation of the <N_i>.
!----------------------------------------
! 2015.12.14
! In the calculation of Kapp, we must recall that the concentration of the 
! current cluster size is used in tandem with the concentration of the monomers
! such that: Kapp(i) = conc(i)/conc(1)**i
!
! There are two ways that we concidered formulating the weighting parameter
! for each Kapp in the geometric mean:
!       ---------------------------------
!       1) Use the relative weight of the standard deviation in the cluster
!       frequency such that each entry to the geometric mean calculated from
!       Kapp is weighted relative to the data quality of that cluster sizes.
!       The weighting is thus the standard deviation in the cluster frequency
!       divided by said cluster frequency for cluster size i: 
!       NrelW(i)=Nstdv(i)/Navg(i)       (w(i)=1/NrelW(i))
!       REVISION: now NrelW(i)=Navg(i)/Nstdv(i) (w(i)=NrelW(i))
!       ---------------------------------
!       2) Use the relative weight based on the two standard deviations
!       (relative to their respective cluster frequencies) such that for a
!       cluster size i:
!       KrelW(i)=NrelW(i)/(NrelW(1)**i) (w(i)=1/KrelW(i))
!       ---------------------------------
! The idea behind the weighting scheme is to weight the contribution from a 
! data point to the geometric mean for a particular cluster size based on the 
! quality of that data point. 
! In method (2), due to raising the relative
! weighting of the monomer's contribution to the cluster size, this weighting
! tends to over emphasize the monomer's relative standard deviation and as such
! leads to weighting the data sets in their entirity based on the quality of a 
! single data point (the monomer) rather than the data point in question. The
! geometric mean from said weighting shows an unequivocal preference for a
! single data set when tested on MTBE simulations with cluster sizes from 
! i=1-30 and then a sharp drop down to the remaining data
! sets for i>30.
! Method (1) works better in the sense that each data point has an individual
! weight (but does not account for the monomer concentrations data quality in
! the generation of Kapp) and returns a curve that is rather similar to the
! unweighted geometric mean.
! 
write (*,*) "##########"
write (*,*) "# GENERATING THE WEIGHTED GEOMETRIC MEAN OF EQUILIBRIUM ASSOCIATION CONSTANTS, K"
write (*,*) "##########"

open(UNIT=10,FILE="MS-multiset-xKguess-GeoMeanWeighted.xvg",form='FORMATTED',status='REPLACE',action='WRITE')
write (10,*) startTime
write (10,*) "## XVG file of the initial guesses for Kactual based on Kapparent"
write (10,*) "## using a weighted geometric mean where the weighting is "
write (10,*) "## NrelW(i) the inverse of the the relative standard deviation in"
write (10,*) "## the frequency that an i-mer occurs in simulation."
write (10,*) "##        NrelW(i)=aveN(i)/stdvN(i)"
write (10,*) "##        Kactual(i)=(Sum_iset^nsamp log(Kapp(i,iset))*"
write (10,*) "##                   NrelW(i,iset))/"
write (10,*) "##                   (Sum_iset^nsamp NrelW(i,iset)) "
write (10,*) "## COL | DESCRIPTION" 
write (10,*) "##  1  | nsamp, number of data points averaged over."
write (10,*) "##  2  | i, the cluster size."
write (10,*) "##  3  | Kactual(i), the equilibrium rate constant for an i-mer."
write (10,*) "##  4  | Delta_G(i), the change in gibbs free energy for the "
write (10,*) "##       i-mer with Kactual(i) and a monomer concentration of "
write (10,*) "##       of ",nc1," molecules/nm^3."
write (10,*) "## "

Kactual(1)=1.0d0
nsamp=0
do i=2,Nmax  
  if (Kactual(i).eq.0.0d0) then  ! that is, if it hasn't been pre-set
    do iset=1,Nsets
      if (ci_gt_0(i,iset)) then
        nsamp(i)=nsamp(i)+1
      endif
    enddo
    if (weightSum(i).ne.0.0d0) then
      Kactual(i) = 0.0d0
      do iset=1,Nsets
        if (ci_gt_0(i,iset)) then
          Kactual(i)=Kactual(i)+NrelW(i,iset)*log(Kapp(i,iset))
        endif
      enddo
      Kactual(i)=exp(Kactual(i)/weightSum(i))
      write (10,FMT3) nsamp(i), i, Kactual(i), -log(Kactual(i))-(i-1)*log(nc1)
    else  
      Kactual(i)=0.0d0
    endif
  endif
enddo
close(unit=10)

write (*,*) "##########"
write (*,*) "# STARTING ITERATIONS OF FITTING"
write (*,*) "##########" 

call fittingRoutine (Nmax, NsetMax, startTime, niter, alpha, Nsets, aveN, &
aveNout, aveNconstrained, nsamp, Nsim, NsimMax, nci_gt_0, ci_gt_0, NrelW, &
weightSum, V, Kapp, Kactual, ilist, revlist, most_prob_comp, nc1)


write (*,*) "##########"
write (*,*) "# FINISHED FITTING"
write (*,*) "##########"

open(unit=3,file='MS-multiset-xOUTPUT.xvg',form='FORMATTED',status='REPLACE',action='WRITE')
write (3,*) startTime
write (3,*) "# Multiset output"
write (3,*) "# COL | VARIABLE"
write (3,*) "#  1  | Cluster size, i"
write (3,*) "#  2  | Concentration via fit (N_i/V) in molec/nm^3"
write (3,*) "#  3  | Concentration via input (N_i/V)"
write (3,*) "#  4  | K_actual: Keff from the fitting"
write (3,*) "#  5  | K_apparent: Keff from the input"
write (3,*) "#  6  | N_i at most probable composition."
write (3,*) "#  7  | Delta G based on K_actual and a monomer concentration of"
write (3,*) "#     | ",nc1," molec/nm^3"
write (3,*) "#  "
do iset=1,nsets
  write (iterString,"(I0.5)") iset
  open(unit=4,file='MS-multiset-xDataSet'//iterString//'-OUTPUT.xvg',&
  form='FORMATTED',status='REPLACE',action='WRITE')
  write (4,*) "# Multiset output for data set ", iset
  write (4,*) "# COL | VARIABLE"
  write (4,*) "#  1  | Cluster size, i"
  write (4,*) "#  2  | Concentration via fit (N_i/V) in molec/nm^3"
  write (4,*) "#  3  | Concentration via constrained fit (N_i/V) in molec/nm^3"
  write (4,*) "#  4  | Concentration via input (N_i/V)"
  write (4,*) "#  5  | K_actual: Keff from the fitting"
  write (4,*) "#  6  | K_apparent: Keff from the input"
  write (4,*) "#  7  | Delta G based on K_actual and a monomer concentration of"
  write (4,*) "#     | ",nc1," molec/nm^3"
  write (4,*) "#  8  | Data's weighting: NrelW(i,iset)/weightSum(i)"
  write (4,*) "#  "
  write (3,*) "#-------------------------"
  do i=1,Nsim(iset)
    if (ci_gt_0(i,iset)) then
      write (3,FMT2) i,aveNout(i,iset)/V(iset),&
      aveN(i,iset)/V(iset),Kactual(i),&
      Kapp(i,iset),most_prob_comp(i,iset), -log(Kactual(i))-(i-1)*log(nc1)
      write (4,FMT5)i,aveNout(i,iset)/V(iset),aveNconstrained(i,iset)/V(iset),&
      aveN(i,iset)/V(iset),Kactual(i),&
      Kapp(i,iset),-log(Kactual(i))-(i-1)*log(nc1),&
      NrelW(i,iset)/weightSum(i)
    endif
  enddo
  close(UNIT=4)
enddo
close(UNIT=3)

open(UNIT=13,FILE="MS-multiset-xKactual.xvg",form='FORMATTED',status='REPLACE',action='WRITE')
write (13,*) startTime
write (13,*) "# This file contains the association constants from the fit"
write (13,*) "# in a format that is designed for implementing a continuation"
write (13,*) "# in the fitting with these values as the initial guesses."
write (13,*) "# "
write (13,*) "# COL | VARIABLES"
write (13,*) "#  1  | Cluster size"
write (13,*) "#  2  | Dummy variable 1"
write (13,*) "#  3  | Dummy variable 2"
write (13,*) "#  4  | Association constant after fitting"
write (13,*) "#  5  | Corresponding Delta_G for a monomer concentration of "
write (13,*) "#     | ",nc1," molec/nm^3: -log(Kactual(i))-(i-1)*log(",nc1,")"
write (13,*) "# "
do i=1,NsimMax
  write (13,*) i,Kactual(i),-log(Kactual(i))-(i-1)*log(nc1)
enddo
close(UNIT=13)




!----------------------------------------
! GETTING THE DATE AND TIME
!----------------------------------------
call date_and_time(DATE=date,TIME=time,ZONE=zone)
finishTime="# Date: "//date(1:4)//"."//date(5:6)//"."//date(7:8)//&
"\n # Time: "//time(1:2)//":"//time(3:4)//":"//time(5:10)//&
" (UTC"//zone(1:3)//":"//zone(4:5)//")\n"
call cpu_time(finishTime_cpu)

write(*,*) "START time: ",startTime
write(*,*) "FINISH time: ",finishTime
write(*,*) "Time elapsed = ",(finishTime_cpu-startTime_cpu)/(60.0*60.0)," hrs"

end  ! end main

!----------------------------------------
!--------------SUBROUTINES---------------
!----------------------------------------

real*8 function Pcalc(Psave,facinv,q,comp,Nsim,Nmax,prestart)
  implicit none
  real*8:: facinv(0:Nmax),q(0:Nmax),Psave(0:Nmax+1)
  integer:: comp(0:Nmax+1,3)
  integer:: Nsim, ipoint,Nmax,prestart
  Pcalc=Psave(prestart)
  ipoint=comp(prestart,3)  ! working down from largest unchanged aggregate 
  do while (ipoint.gt.0)
    Psave(ipoint)=Pcalc*(q(ipoint)**(comp(ipoint,1)))*facinv(comp(ipoint,1))
    Pcalc=Psave(ipoint)
    ipoint=comp(ipoint,3)
  enddo
  Psave(0)=Pcalc
  if (Pcalc.lt.0.0d0) then
    Pcalc=1.0d0
    ipoint=comp(Nsim+1,3)
    do while (ipoint.gt.0)
      write (*,*) ipoint,q(ipoint),comp(ipoint,1),Pcalc,facinv(comp(ipoint,1))
      Pcalc=Pcalc*(q(ipoint)**(comp(ipoint,1)))*facinv(comp(ipoint,1))
      ipoint=comp(ipoint,3)
    enddo
    stop
  endif
return
end

subroutine factorial_init (facinv,Nmax)!(factorial, facinv, Nmax)
  implicit none
  real*8:: factorial(0:Nmax),facinv(0:Nmax)
  integer:: Nmax,i
  ! note: floating-point error for Nmax > 170
  ! 2015.12.11 Experimentation using the Exp(Log()) trick has shown that the
  ! sum of log(i)'s remains within the floating-point range of sizes. Once the
  ! exponent of that sum is taken however, we incur the same floating-point 
  ! error for Nmax > 170. This seems non-negotiable unless there's a way to do
  ! all the rest of the calculations using the Exp(log()) format and only take 
  ! the exponent at the very end. 
  factorial(0)=1.0d0
  facinv(0)=1.0d0
  do i=1,Nmax
    factorial(i)=factorial(i-1)*i
    facinv(i)=1.0d0/factorial(i)
  enddo 
  return
end

subroutine forward (q,facinv,Nsim,ilist,revlist,nci_gt_0,Nmax,&
aveNout,most_prob_comp,ncomp)
!(q,facinv,aveNout,Nsim,Nmax,ilist,revlist,nci_gt_0,most_prob_comp,ncomp)
  ! input   : q, facinv, Nsim, ilist, revlist, nci_gt_0, Nmax
  ! returns : aveNout, most_prob_comp, ncomp
  implicit none
  real*8:: q(0:Nmax),facinv(0:Nmax),aveNout(0:Nmax)
  real*8:: Qtot,P,Psave(0:Nmax+1)
  real*8:: Pcalc,biggest_P
  integer:: Nmax,Nsim,nci_gt_0
  integer:: ilist(0:Nmax+1),revlist(0:Nmax+1)
  integer:: inext,irest,ibig,nsplit
  integer:: ipoint,iloop,inext_minus_1,irestnext,inextpoint,prestart
  integer:: comp(0:Nmax+1,3)  ! comps with linked list attached
  !         comp(i,1) contains the number of i-mers in the composition/partition
  !         comp(i,2) contains the next-highest value of j with comp(j,1)>0
  !         comp(i,3) contains the next-lowest value of j with comp(j,1) >0
  !         comp(Nsim+1,3) = the largest aggregate present in the composition
  !         comp(0,2) = the smallest aggregate (incl. free monomers) present
  !         comp(i,3) =0 if i is the smallest aggregate
  !         comp(i,2) =Nsim+1 if i is the greatest aggregate
  integer:: most_prob_comp(0:Nmax+1),ncomp
  ! generate all possible compositions of a system with Nsim monomers
  ! start off with composition containing maximum number of biggest 
  ! defined clusters
  
  ibig=ilist(nci_gt_0)
  inext=ibig
  
  aveNout=0.0d0
  comp=0
  comp(inext,1) = Nsim/inext ! has to do with the sorting of partitions. 
  ! want to find the largest integer number of ibig in the system as the first
  ! partition.
  comp(inext,2) = Nsim+1  ! linked list pointing upwards; 
                          ! index corresponds to ilist
                          ! Defines the end of the list so to speak.
  comp(Nsim+1,3)=inext
  irest=mod(Nsim,inext)
  do while (irest.ne.0) !! JUST AN INITIATION OF THE HIGHEST OCC i-mer PARTIT.OR
    irestnext=ilist(revlist(irest))  ! skip down in case irest-mers are not present
    comp(irestnext,1)=irest/irestnext
    comp(irestnext,2)=inext
    comp(inext,3)=irestnext
    !write (*,*) inext,irest,irestnext
    irest=mod(irest,irestnext)
    inext=irestnext
  enddo   
  
  comp(0,2)=inext !Points up in partition
  comp(inext,3)=0 !Points down in partition
  inext=max( comp(1,2),comp(0,2)) !Lowest non 1-mer i-mer that is occupied. 
  
  Qtot=0.0d0
  Psave=1.0d0
  
  comp(0,1)=0
  
  prestart=Nsim+1  
  
  ncomp=0
  biggest_P=0.0d0
  do while (inext.le.Nsim)
    ! calculate new probability
    ! weight contributions from i-mers with i >= prestart are
    ! unchanged, save time by not recalculating them
    P=Pcalc(Psave,facinv,q,comp,Nsim,Nmax,prestart)
    if (P.gt.biggest_P) then
      biggest_P=P
      most_prob_comp=comp(:,1)
    endif
    
    ! count current composition
    ncomp=ncomp+1
    Qtot=Qtot+P
    ipoint=comp(Nsim+1,3)
    do while (ipoint.gt.0)
      aveNout(ipoint)=aveNout(ipoint)+comp(ipoint,1)*P
      ipoint=comp(ipoint,3)
    enddo
    
    ! determine lowest aggregate size greater than a monomer
    inext=max( comp(1,2),comp(0,2))
    ! explanation: either 0 points to 1 and 1 points up to the next-highest 
    ! size, or there are no free monomers and 0 points directly to inext.
    ! divide up all the monomers belonging to one inext-mer, with free monomers added on, ! into the greatest possible aggregate size (i.e., ilist(revlist(inextpoint-1))- mers )

    prestart=comp(inext,2)
    inextpoint=revlist(inext)
    inext_minus_1=ilist(inextpoint-1) ! note that this will be inext-1 or, if conc(inext-1)=0, the next lower aggregate size
    nsplit=(inext+comp(1,1))/inext_minus_1
    irest=mod(inext+comp(1,1),inext_minus_1)
    comp(1,:)=0
    comp(0,:)=0
    comp(inext,1)=comp(inext,1)-1
    comp(inext_minus_1,1)=nsplit
    if (comp(inext,1).gt.0) then  ! link up and down to inext_minus_1
      comp(inext,3)=inext_minus_1
      comp(inext_minus_1,2)=inext
    else                          ! if inext is cleared out, shift links down
      comp(comp(inext,2),3)=inext_minus_1
      comp(inext_minus_1,2)=comp(inext,2)
      comp(inext,2)=0; comp(inext,3)=0
    endif
    
    ! deal with remainder
    iloop=inext_minus_1
    do while (irest.ne.0)
      irestnext=ilist(revlist(irest))  ! skip down in case irest-mers are not present
      comp(irestnext,1)=irest/irestnext
      comp(irestnext,2)=iloop
      comp(iloop,3)=irestnext
      irest=mod(irest,irestnext)
      iloop=irestnext
    enddo
    comp(iloop,3)=0
    comp(0,2)=iloop
  enddo  ! finished one set of compositions for one trial set of q(i)
  aveNout=aveNout/Qtot
  return
end

subroutine fittingRoutine (Nmax, NsetMax, startTime, niter, alpha, Nsets, &
aveN, aveNout, aveNconstrained, nsamp, Nsim, NsimMax, nci_gt_0, ci_gt_0, &
NrelW, weightSum, V, Kapp, Kactual, ilist, revlist, most_prob_comp, nc1)
  implicit none
  
  character(100) :: startTime
  character(5) :: iterString ! String representation of an integer for file name

  integer :: FN_CONV, FN_ITER !! Unit numbers for writing files
  integer :: PF ! Printing frequency
  integer :: Nmax ! Maximum i-mer size
  integer :: NsetMax ! Maximum number of data sets
  integer :: niter ! Number of iterations
  integer :: Nsets ! Number of data sets
  integer :: iter, iset, i ! Counters
  integer :: nsamp(0:Nmax) ! Number of samples for each i-mer
  integer :: Nsim(Nsetmax), NsimMax 
  integer :: ncomp(Nsetmax)
  !! its part of the forward subroutine but isn't used in the main program. Not
  !! sure its critical that it comes back out of the subroutine.
  integer :: nci_gt_0(Nsetmax)
  integer :: ilist(0:Nmax+1,Nsetmax), revlist(0:Nmax+1,Nsetmax)
  integer :: most_prob_comp(0:Nmax+1,Nsetmax)

  real*8 :: aveN(Nmax,Nsetmax), aveNout(0:Nmax,Nsetmax)
  real*8 :: aveNconstrained(0:Nmax,Nsetmax)
  real*8 :: NrelW(Nmax,Nsetmax), weightSum(Nmax)
  real*8 :: V(Nsetmax)
  real*8 :: Kapp(Nmax,Nsetmax), Kactual(Nmax)
  
  real*8 :: alpha, Kapp_current, nc1
  real*8 :: facinv(0:Nmax), q(0:Nmax)
  real*8 :: Kfac, aveNoutC(0:Nmax)
  real*8 :: converCrit(Nmax),converCritTot,converVar
  real*8 :: converCritGM(Nmax),converCritTotGM

  integer :: nci_gt_0C, ncompC, NsimC
  integer :: ilistC(0:Nmax+1), revlistC(0:Nmax+1), most_prob_compC(0:Nmax+1)
  
  logical :: ci_gt_0(Nmax,Nsetmax)
  ! Bools that control whether or not to print in that iteration.
  ! Evaluated only once at the begining of the iteration.
  logical :: printing, printing_short
 
  FN_CONV=12
  FN_ITER=14
  aveNoutC=0.0d0

  ! Defining the factorials and inverse factorials using the module/function
  ! factorial_init
  call factorial_init (facinv,Nmax)

  open(UNIT=FN_CONV,FILE="MS-multiset-xConvergenceCriteria.xvg",&
  form='FORMATTED', status='REPLACE',action='WRITE')
  write(FN_CONV,*) startTime
  write(FN_CONV,*) "# This file contains a plottable progression for the "
  write(FN_CONV,*) "# convergence criteria, primarily the result of sum: "
  write(FN_CONV,*) "#   ((n_model(i,iset)-n(i,iset))**2/n_stdv(i,iset)**2)"
  write(FN_CONV,*) "# over all the data. The desired trend should be to see"
  write(FN_CONV,*) "# the value of this go down and the iteration step goes"
  write(FN_CONV,*) "# up."
  write(FN_CONV,*) "# COL | VARIABLE "
  write(FN_CONV,*) "#  1  | iteration step"
  write(FN_CONV,*) "#  2  | sum_(i) ( prod_(iset) ("
  write(FN_CONV,*) "#     | ((n_model(i,iset)-n(i,iset))**2/n_stdv"
  write(FN_CONV,*) "#     | (i,iset)**2)**(NrelW(i,iset)/weightSum(i))"
  write(FN_CONV,*) "#     | ) )"
  write(FN_CONV,*) "# "

  do iter=1,niter
    ! Defining the printing frequency.
    printing=((mod(iter,PF).eq.0).or.(iter.eq.1))
    printing_short=.TRUE.
    if (printing) then
      write(*,*) "Iteration ",iter
    endif
    if (iter.eq.1) then
      ! Redefining the printing frequency for the post constrainment fitting.
      PF=niter/20
      printing=(mod(iter,PF).eq.0).or.(iter.eq.1)
      write(*,*) "#-----------------------\n # niter = ",niter,&
      "\n # alpha = ",alpha,&
      "\n # The frequency with which data is printed will be every",&
      PF,"\n # steps. This amounts to a total of 20 prints over the ",&
      niter,"\n # total iterations."
    endif
    do iset=1,Nsets
      NsimC=Nsim(iset) ! total number of monomers in this run
      ilistC=ilist(:,iset)
      revlistC=revlist(:,iset)
      nci_gt_0C= nci_gt_0(iset)
      q(0)=1.0d0
      q(1)=1.0d0
      do i=1,Nmax
        ! generate "q" (actually scaled to V) array from current best
        ! guess Kactual
        !----------------------------------------
        ! INITIAL GUESS OF q(i) FROM APPARENT K
        ! First guess at finding single-aggregate partition functions relative
        ! to monomer. (For a single data set, this will set Kactual = Kapp)
        !----------------------------------------
        ! Note: q0(i) is the standard state partition function for an i-mer.
        !       V0 is the standard state volume, assumed to be equal to 1.
        !       K(i) is the equilibrium constant for the i-mer and equals
        !               q0(i)/q0(1)**i
        ! q(i)=q0(i)*(V/V0)
        !     =K(i)*(q0(1)**i/q0(i))*q0(i)*(V/V0)  
        !               !multiplying by 1=K(i)*K(i)**-1
        !     =K(i)*q0(1)**i*(V/V0)
        !     =K(i)*q(1)**i*(V0/V)**i*(V/V0)       !q0(1)=q(1)*(V0/V)
        !     =K(i)*q(1)**i*(V0/V)**(i-1)
        !     =K(i)*q(1)**i*V**(1-i)               !V0=1
        !     =K(i)*V**(1-i)                       !q(1)=1.0d0
        !----------------------------------------
        q(i)=Kactual(i)*V(iset)**(1-i)
      enddo
      ! Subroutine forward returns aveNoutC, most_prob_compC, and ncompC
      call forward (q,facinv,NsimC,ilistC,revlistC,nci_gt_0C,Nmax,&
      aveNoutC,most_prob_compC,ncompC)
      ! Generate new set of <n_s> from current Kactual list
      aveNout(:,iset)=aveNoutC
      most_prob_comp(:,iset)=most_prob_compC
      ncomp(iset)=ncompC
    enddo
    ! Tracking the convergence of Kactual
    ! The following calculation keeps track (in a loose sense) of the
    ! quality of the fitting for each i-mer individually and also the
    ! data set as a whole. Contributions here are again weighted by
    ! NrelW/weightSum.
    converCritGM=0.0d0
    converCritTotGM=0.0d0
    do i=2,nmax
      do iset=1,nsets
        if (ci_gt_0(i,iset)) then
          converVar=(NrelW(i,iset)/weightSum(i))*&
          log((aveNout(i,iset)-aveN(i,iset))**2/aveN(i,iset)**2)
          converCritGM(i)=converCritGM(i)+converVar
        endif
      enddo
      if (converCritGM(i).ne.0.0d0) then
        converCritGM(i)=exp(converCritGM(i))
      endif
      converCritTotGM=ConverCritTotGM+converCritGM(i)
    enddo
    if (printing_short) then ! calculate total convergence and print
      converCritTotGM=converCritTotGM/dfloat(NsimMax)
      write(FN_CONV,*) iter, converCritTotGM
    endif
    !! END OF CONVERSION CRITERIA CALCULATION
    do i=2,nmax
      Kfac=0.0d0
      if (nsamp(i).gt.0) then
        do iset=1,nsets
          if (ci_gt_0(i,iset)) then
            ! calculate the apparent association constant from the latest
            ! iteration. Make a new guess based on the ratio of the current 
            ! constant to the true value; each set was originally weighted 
            ! equally:
            !     Kfac=Kfac*(Kapp(i,iset)/Kapp_current)**(alpha/dfloat(nsamp))
            ! but that has been changed using NrelW/weightSum which
            ! should appropriate a greater weight to data points that have a
            ! smaller relative standard deviation.
            !Kfac=Kfac+log(aveN(i,iset)/aveNout(i,iset))*(alpha*&
            !NrelW(i,iset)/weightSum(i))
            Kapp_current=(V(iset)**(i-1))*aveNout(i,iset)/(aveNout(1,iset)**i)
            Kfac=Kfac+log(Kapp(i,iset)/Kapp_current)*(alpha*&
            NrelW(i,iset)/weightSum(i))
          endif
        enddo
        Kactual(i)=Kactual(i)*exp(Kfac) 
      endif
    enddo
    Kactual=Kactual/Kactual(1) ! Normalized with respect to Kactual(1)=1
    if (printing) then
      write(iterString,"(I0.5)") iter/PF
      open(UNIT=FN_ITER,FILE="MS-multiset-xIter"//iterString//"Kactual.xvg",&
      form='FORMATTED',status='REPLACE',action='WRITE')
      write(FN_ITER,*) startTime
      write(FN_ITER,*) "# Iteration ",iter
      write(FN_ITER,*) "# i, Kactual(i), DeltaG(i)"
      write(FN_ITER,*) "# "
      do i=1,NsimMax
        write (FN_ITER,*) i, Kactual(i), -log(Kactual(i))-(i-1)*log(nc1)
      enddo
      close(UNIT=FN_ITER)
    endif
  enddo
  close(UNIT=FN_CONV)
  return
end
