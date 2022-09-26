 --- renumber residues of pdb file from pscdb pdb 
 --- to make a set of ligand bind pdb
 --- estimate alignment of residues from xyz coordinates
  --- devversion by 2022 Sept Yutaka Ueno AIST Japan

  --- a sample code
  --- bvh-bone data generator 

  
dofile("fit3p.lua")
dofile("pdbmol.lua")

  --- some compat stuff
function getn(t) if(t) then return #t end end


ConvRes2Atom={}

function main()
	print(arg[0],arg[1])
---if(1) then return end

									--- 
	PDBNAME="1bbu"
	RFILE="p048150.pdb" --- renumbered pdb (with ligand)
	DBCODE="048_1_" --- LYSYL-TRNA SYNTHETASE
	
		--- the data id in PSCDB database
									--- 

  DATADIR="pscdb/"  			--- DATADIR="..\\..\\Downloads\\"

  PCODE="pb0481m"
	  --- but this time take another one
	
	BFILE="p8-"..PDBNAME..".bvh"

	io.output(BFILE)

		--- reference pdb 
	io.input(RFILE)
	refpdb=msReadPdb()
	if(not refpdb) then
		print("reference pdb error",RFILE)
		return
	end
	print("reference pdb",RFILE,#refpdb)

	--- take only the first chain including hetatom
	--- please check if the target hetatom is in an another chain
	cname=refpdb[1].cname
	for ia=#refpdb,1,-1 do
		if(refpdb[ia].cname~=cname) then --- or ( refpdb[ia].het) then 
			--table.remove(refpdb,ia)
		end
	end

	print("first chain and hetatom",cname,#refpdb)

		--- load moprhing pdb data file
			 --- support full atom file and only backbone files
			 --- only the first chain
	selectchain=nil
	tfs={49,40,30,20,10,1} --- ,30,20,10,0}    ---{1,10,20,30,40,49}  --- tiime frame number
	--tfs={1,2,3,4}
			 --- load mutilple files
	allpdb={}
	ligand={}
	for itm=1,getn(tfs) do
		hetatmlist=""
		hetatm=""
		if(tfs[itm]>=0) then
			TFILE=DATADIR..DBCODE..string.format("%0.3d",tfs[itm])..".pdb"
			--TFILE=string.format("%d",tfs[itm])..".pdb"
			OUTFILE=PCODE..string.format("%0.3d",tfs[itm])..".pdb"
		else
			TFILE=DATADIR..DBCODE..".pdb"
		end
		print(itm,TFILE)
		io.input(TFILE)

		selected={}
		mol=msReadPdb()

		selectchain=mol[1].cname


		--- get all atoms
		for ia=1,getn(mol) do
			atom=mol[ia]
				--- there are only alpha carbons
			if(renumbered) then
				atomref=renumbered[ia]
			else
				atomref,dist=findinpdb(refpdb,atom)
				if(dist>1) then
					print("far atom match",atom.resname,atom.resno, atom.aname,atomref.resname,atomref.resno)
				end
			end
			atom.resname=atomref.resname
			atom.resno=atomref.resno
			atom.cname=atomref.cname

		end

		  --- get additional ligand model with HETATM
		  --- static model
		  --- if(1) then return end ------
		print("save",OUTFILE)
		io.output(OUTFILE)
		msWritePdb(mol)
		--- use previous pdb for the ref pdb
		if(not renumbered) then
			renumbered=mol
		end

	end
end

  --- get clothest atom in the refpdb
function findinpdb(refpdb,atom)
	local atomref
	local mindist

	for ia=1,#refpdb do
		btom=refpdb[ia]
		dx=btom.x-atom.x
		dy=btom.y-atom.y
		dz=btom.z-atom.z
		dist=dx*dx+dy*dy+dz*dz
		if (not mindist) or (mindist>dist) then
			if(atom.aname==btom.aname) then
				mindist=dist
				atomref=btom
			end
		end
	end
	return atomref,math.sqrt(mindist)
end


--------------------



function debug3pt(firstp,midp,lastp,firstp2,midp2,lastp2,mop)
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	  midp.x,midp.y,midp.z,midp2.x,midp2.y,midp2.z))
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	  firstp.x,firstp.y,firstp.z,firstp2.x,firstp2.y,firstp2.z))
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	  lastp.x,lastp.y,lastp.z,lastp2.x,lastp2.y,lastp2.z))
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	mop.mtx[1], mop.mtx[2], mop.mtx[3], mop.mtx[4], mop.mtx[5], mop.mtx[6]))
	print(string.format("%8.3f %8.3f %8.3f",
	mop.mtx[7], mop.mtx[8], mop.mtx[9] ))
end

--------------------
main()
 
