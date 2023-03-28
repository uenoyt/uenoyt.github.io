MOSBY : Molecular Structure Browser 
ASHLEY: Application Support Hybrid Library for EasY programming

developer release  0.936			Jan 15, 2003
						Yutaka Ueno,
						Electrotechnical Laboratory
------------------
About This Release
------------------

MOSBY is a molecular structure viewer program for the use in structural 
biology and protein structure prediction studies. ASHLEY is a programing 
library for the plug-in software configuration. This distribution of the 
source code is to demonstrate a new plug-in software architecture that is 
described in a published articles:

 Yutaka Ueno and Kiyoshi Asai. "A New Plug-in Software Architecture 
 Applied for a Portable Molecular Structure Browser", in the Proceeding, 
 Intelligent System for Molecular Biology 1997(Gaasterland Eds.) 329-332.

This is a beta version tested by various users. We have delayed 
our official release until users community required to fix.
We are glad to provide this release and evaluated by anybody interested
in our programs and its new software method.
Any comment would be appreciated to ueno@etl.go.jp.
Please refer the adobe paper in works with the use of this software.

-------------------
MOSBY
-------------------

makeup/		a final binary distribution
main/		mosby source code 

	mosby.c		main command processing
	molexp.c	main molecular list management
	mosby.ui	GUI definitions
	demo.pdb	Adenosin 3 phosphate

	club2.c		"club" application data list support
	dbox.c		tipical dialog box for application
	edge.c		edge line drawing module
	geomtx.c	more 3x3 matrix for geometry calc
	gotolua.c	very bad LUA language file reader 

	moltype.h
	molmem.c	molecular data list
	molnum.c	counting up moleucles

	molgrp.c	molecular group 
	molsel.c	select molecule
	molinf.c	chemical information
	molpdb.c	pdb file format filter module
	mosfile.c	mosby molecular experiment file filter module
	molmdt.c	modeling translation for symmeotric unit
	molxtal.c	crystallographic symmetry operation
	molview.c	display moleucle
	molvis.c	molecular representaion 
	mosdsp.c	display list for the molecular representation
	mosedit.c	modify atomic data
	mosxyz.c	general coordinate array
	mouse3d.c	mouse manipulation 
	mvsqa.c		structure table window
	pluginfo.c	plugin information dialog
	qsphere.c	copy template sphere by bitblt with depth
	vw3d.c		3d view window
	yfilter.c	file format dispatcher
	yftn.c		fortran format file read utility


-------------------
PLUG-INS 
-------------------
edmap/		electron density map
filters/	MDL file filter
fitmol/		docking study plug-in template
fsample/	FORTRAN plugin sample
peach/		molecular dynamics trajectory
gnuplot/	gnu-plot plugin experimental module
rasmol/		let's bind part of RasMol code
sample/		sample plug-in 
tagl/		texture mapping experimental module
tmv/		make a TMV helix
str4sgi/	sgi stereo hardware enabler

-------------------
NOTE 
-------------------


---------------
CONTACT
---------------
Yutaka Ueno
Computational Biology Research Center
AIST Tokyo
uenoyt@ni.aist.go.jp


/*///////////////////////////////////////////////////////////////////
Copyright (c) 2001 National Institute of Advanced Industrial 
		   Science and Technology (AIST)

Permission to use this material for non comercial purpose is hereby
granted, provided that the above copyright notice and this permission
notice appear in all copies. A registration is required to modify,
copy, and distribute this material.
AIST AND AUTHOR MAKE NO REPRESENTATIONS ABOUT THE ACCURACY OR
SUITABILITY OF THIS MATERIAL FOR ANY PURPOSE. IT IS PROVIDED "AS IS",
WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
///////////////////////////////////////////////////////////////////*/

