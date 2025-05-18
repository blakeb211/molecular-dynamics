! Copyright (C) 2025 Blake Baird
! See https://factorcode.org/license.txt for BSD license.
! Molecular Dynamics with a pairwise potential function 
! Visualize it in some way

USING: math random sequences grouping namespaces 
       math.vectors kernel splitting columns math.constants
       arrays math.functions compiler.utilities prettyprint
       vocabs.loader ;
IN: md

SYMBOL: BOXLEN
SYMBOL: TAPE
SYMBOL: DT

! posx posy posz velx vely velz
: init-box ( n -- seq )
[ 
  3 [ BOXLEN get random ] replicate 
  { 0.0 0.0 0.0 } append
] replicate ;

: save-state ( x -- x )
dup TAPE get swap suffix TAPE set ;

: init ( n boxlen dt -- x ) 
 dup 1.0 <= [ ] [ assert ] if
 DT set
 BOXLEN set
 { } TAPE set 
 init-box 
 2 [ save-state ] times ;

:: lennard-jones ( dists -- force ) 
 2.0 dists n/v 6.0 v^n :> distsPow6
 dists 2.0 v^n :> distsPow2
 48.0 4.0 * distsPow6 n*v distsPow6 0.5 v-n v* distsPow2 v/ ;

:: calc-force ( index otherpos 6vec -- a(t+dt) )
 otherpos [ [ 3 head 6vec 3 head distance ] 
 [ 3 head 6vec 3 head v- dup norm v/n ] bi 2array ] map 
 dup 0 <column> >array :> dists
     1 <column> >array [ -1.0 v*n ] map :> dirs
     dists lennard-jones 
     ! dup "force" . .
     ! dirs "dirs" . .
     dirs [ n*v ] 2map
     { 0. 0. 0. } [ v+ ] reduce "force" . dup . ; inline

:: velocity-verlet ( 6vec index -- new )
 ! position update
 index TAPE get penultimate nth :> prev6
 6vec 3 head :> r(t)
 ! "prev6" . prev6 .
 ! "6vec" . 6vec .
 DT get :> dt 
 6vec prev6 v- 3 tail :> a(t)
 6vec 3 tail :> v(t)
 ! "a(t)" . a(t) .
 r(t) v(t) dt v*n v+ a(t) 0.5 v*n dt v*n dt v*n v+ :> r(t+dt) 
 ! "r(t+dt)" . r(t+dt) .
 ! force update
 index
 index TAPE get last remove-nth 
 r(t+dt) 
 calc-force :> a(t+dt) 
 ! velocity update
 ! a(t+dt) a(t) v(t)
 a(t) a(t+dt) v+ dt v*n 0.5 v*n v(t) v+ :> v(t+dt)
 r(t+dt) v(t+dt) append ; inline


! iterate velocity and position 
! save state
: next ( state -- state2 ) 
[ velocity-verlet ] map-index ! over each atom
save-state ;

! init system to a well defined 2 atom state
: test-setup ( state -- testState ) 
 drop
 { { { 0.0 0.0 1.5 0.0 0.0 0.0 } { 0.0 0.0 -1.5 0.0 0.0 0.0 } } 
 { { 0.0 0.0 1.5 0.0 0.0 -100.3 } { 0.0 0.0 -1.5 0.0 0.0 +100.3 } } }
 TAPE set
 TAPE get last ;

: iterate ( -- x y z )
"md" reload
 2 5.0 0.0001 init
 test-setup 
 10000 [ next ] times 
 TAPE get
 dup [ 0 swap nth 2 swap nth ] map ;

! xs
! [2.0, 2.25, 8.0]
!
! list(map(LJforce,xs))
! [24.0, -0.1258995652831338, -0.00036603212356567383]
!
! list(map(LJenergy,xs))
! [15.0, 3.399777410804181, -0.00024318695068359375]
