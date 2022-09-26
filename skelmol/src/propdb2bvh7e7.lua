 --- trying to estimate missing motions by adjacent bones
 --- support ligand bound seriese

  --- make a trajectory  motion bones from pdb trajectory 
  ---
  --- the updated coord
  --- tested version by 2022 Yutaka Ueno
  --- AIST Japan

  --- a sample code
  --- bvh-bone data generator 

--- fit to the pdb import add-on in blender
--- x3d yFo zUp
--- bvh yFo zUp 
  
dofile("fit3p.lua")
dofile("pdbmol.lua")

  --- some compat stuff
function getn(t) if(t) then return #t end end


ConvRes2Atom={}

function main()
	print(arg[0],arg[1])
---if(1) then return end

									--- 
		molname="1bbu" 	--- name of pdb for the conformation 1 (ligand free)
		PDBNAME=molname
	DBCODE="pb0481m"	--- the data id in PSCDB database
									--- LYSYL-TRNA SYNTHETASE

---	DBCODE="0000"

  DATADIR="" --"pscdb/"  			--- DATADIR="..\\..\\Downloads\\"


	---  the reference pdb with ligand bind form
	RFILE=DATADIR..molname..".pdb"
	---molname="1gjt1"
	--RFILE="pb0481m050.pdb"

	RFILE=nil
	  --- but this time take another one
	  --- get 2nd st DSSP notation 

	BFILE="pr7e7-"..molname..".bvh"

	io.output(BFILE)

		--- reference pdb 
	if(RFILE) then
		io.input(RFILE)
		refpdb=msReadPdb()
		if(not refpdb) then
			print("reference pdb error",RFILE)
			return
		end
		print("reference pdb",RFILE,#refpdb)


		--- take all chain including hetatom
		--- please check if the target hetatom is in an another chain
		cname=refpdb[1].cname
		for ia=#refpdb,1,-1 do
			if(refpdb[ia].cname~=cname) then --- or ( refpdb[ia].het) then 
				---table.remove(refpdb,ia)
			end
		end
	end

		--- load moprhing pdb data file
			 --- support full atom file and only backbone files
			 --- only the first chain
	selectchain=nil
	tfs={50,49,40,30,20,10,1}    ---{1,10,20,30,40,49}  --- tiime frame number
	--tfs=nil 
	--pdbmdl="mov-1gjta-1gjsa.pdb"
	--tfs={1,2,3,4}
			 --- load mutilple files
	allpdb={}
	ligand={}
	itm=0
	while(true) do
		itm=itm+1
		if(tfs) then
			tfnumb=tfs[itm]
			if( not tfnumb) then
				nframe=itm
				break
			end
		else
			tfnumb=itm
		end

		hetatmlist=""
		if(tfnumb>=0) then
			TFILE=DATADIR..DBCODE..string.format("%0.3d",tfnumb)..".pdb"
		end
		if(pdbmdl) then
			if(itm==1) then
					--- open multi model pdb file
				io.input(pdbmdl)
			end
		else
			print(TFILE)
			io.input(TFILE)
		end

		selected={}
		mol=msReadPdb()

		DEBUG=nil
		if(DEBUG) then
					--- JUST for DEBUG
					for ia=#mol,DEBUG,-1 do
						table.remove(mol,ia)
					end
					for ia=1,#mol do
						print("pd",mol[ia].x,mol[ia].y,mol[ia].aname)
					end
		end

		---print("take the first chain",selectchain, #mol.sequence[selectchain])

		--- get all atoms
		for ia=1,getn(mol) do
			atom=mol[ia]
				--- here only take alpha carbon
			aname=atom.aname
			-- not changing sequence
			cname=atom.cname
			if(cname==" ") then
				cname="A"
				atom.cname=cname
			end
			if(mol[ia].het) then
				if(ligandfix) then
						--- ignore and just take the first locaiton
				else
					table.insert(ligand,atom)
				end
				hetatm=mol[ia].resname
				if(not string.find(hetatmlist,hetatm)) then
					hetatmlist=hetatmlist.." "..hetatm

				end
			else
				table.insert(selected,atom)

				if (aname==" CA ") then
					ConvRes2Atom[atom.resno]=ia
				-- ok
				end
			end
		end

		  --- get additional ligand model with HETATM
		  --- static model

		if(mol) then
			if(#mol==0) then
				--- no more data
				mol=nil
				break
			end

			 --- selected.sequence=mol.sequence
			table.insert(allpdb,mol)
			--- load next MODEL data
			print("pdbone",#mol)
			if(hetatm) then 
				print("hetatom",hetatmlist,#ligand)
				ligandfix=true
			end

			  --- next frame file name 
		else
			mol=nil
			break
		end
	end

	print("hetatom count",#ligand)

		  --- if(1) then return end ------
	if(not refpdb) then
		refpdb=allpdb[1]
	end

	groups=MakeResgroups(refpdb)


function getnamedatom(moldata,cname,resno,aname,start)
	for ia=start or 1,#moldata do
		atomdata=moldata[ia]
		if(atomdata.cname==cname) then
			if(atomdata.resno==resno) then
				if(atomdata.aname==aname) then
					return atomdata
				end
			end
		end
	end
end

--- convert xyz to alpha carbon array at residue index
	--- calc movements vector for each traj
	nframe=getn(allpdb)
	mol=allpdb[1]
	nra=getn(mol)
	print("NRES",nra,"pdb-models",nframe)
	print("chain name in pdb",mol[1].cname)
	print(getn(mol.sequence.A),getn(mol))

	 --- multi chain list
	allchain={}

	for is=1,nframe do
		mol=allpdb[is]
		  -- residue number index xyz[resno]
		ir=0
		ia2=0
		for ia=1,#mol do
					--- print("D",ia,mol[ia].resno,mol[ia].aname)
			resno=mol[ia].resno
			cname=mol[ia].cname
			if(not allchain[cname]) then
				allvec={}
				vec={}
				table.insert(allvec,vec)
				allchain[cname]=allvec
						print("pdb chain",cname)
			end
			allvec=allchain[cname] 
			vec=allvec[is]
			if(not vec) then
				vec={}
				allvec[is]=vec
			end

			resxyz=vec[resno]
			if(not resxyz) then
				resxyz={}
				vec[resno]=resxyz
			end

			atomone=mol[ia]
			if(mol[ia].aname==" CA ") then
				resxyz.caatom=mol[ia]
				resxyz.caidx=ia
				  -- make a residue data
				--- use the sparse array for residue number index
				--- use updated residue number
				resxyz[1]=atomone.x
				resxyz[2]=atomone.y
				resxyz[3]=atomone.z
			elseif(mol[ia].aname==" O  ") then
				resxyz.oatom=mol[ia]
			elseif(mol[ia].aname==" C  ") then
				resxyz.c_atom=mol[ia]
			---debug  if(ia<111) then print("resxyz",ir,resxyz[1],resxyz[2],resxyz[3],ia,mol[ia].resno) end  				
			end
		end
	end

if(1) then
	pdbchain={}
	--- delete frame 1 
	for idx,cnt in pairs(allchain) do
		allvec=cnt
		pdbchain[idx]=cnt[1]
		if(getn(allvec)>1) then
			table.remove(allvec,1)
			print("del frame 1",idx)
		end
	end
end

--------------------same as 

	  --- calc gravity and vector for groups from sst (secondary structure )
	bonelist={}
	for ig=1,getn(groups) do
		span=groups[ig]
		  --- support multi chain sequence
		allvec=allchain[span.cname]
		  --- set 1st, last, and mid-control points
		  --- for helix or sheet, apply long range bone
		  --- the control point smid
		smid=math.floor((span[2]+span[1])/2)

		-- special for the alpha helix, an emprical tip
		---if(sst=="H")or(sst=="E") then
		if(span.sst=="H") then
			smid=span[1]+5
		else
			--- smid=span[1]+1 ---4
		end
----trial
----		if(span.sst=="S") then
----			smid=span[1]+1
----		end
		---
		len=span[2]-span[1]+1
		isp1=span[1]
		isp2=span[2]
		iep1=span[1]
		iep2=span[2]
		if(len==3) then
			isp2=isp1+1
			iep1=iep2-1
		elseif(len<=2) then
			isp2=isp1
			iep1=iep2
		else
			isp2=isp1+2
			iep1=iep2-2
		end
		  --- for secondary structure groups, apply
		  --- this unique rule for the 3-point
		if(len==6) then
			isp2=isp1+1
			iep1=iep2-1
		end

		  --- other loop regions, divide into 2 bodies
		if(len>2) then
			if(span.sst~="H")and(span.sst~="E") then
				isp2=smid-1
				iep1=smid+1
			end
		end

		span.smid=smid
		span.isp1=isp1
		span.isp2=isp2
		span.iep1=iep1
		span.iep2=iep2
	end

	  --- duplicate the groups for dimer no?


	  --- calc average resxyz coordinates for each spans and motion vector

	for ig=1,getn(groups) do
		span=groups[ig]

		  --- print("debug span",ig,span[1],span[2],span.name)
		allvec=allchain[span.cname]
		  -- print("Gr",span.cname,span[1],span[2],span.smid)
		 --- set bone coordinate from atom resxyz
		atomxyz=pdbchain[span.cname]

		  --- for single residue group
		if(span[1]+1==span[2]) then
			if(atomxyz[span[2]]) then
				span.grav=XYZcalcavg(atomxyz,span[1],span[2])
				span.gvec,log =XYZveccalcavg(allvec,span[1],span[2], span.grav)  -- if atom is missing use default: span.grav
			else
				span.grav=XYZcalcavg(atomxyz,span[1],span[1])
				span.gvec,log =XYZveccalcavg(allvec,span[1],span[1], span.grav)
			end
			if(not span.grav) then
					print("warning this pdb file has no atoms for group ",ig,span[1],span[2],span.cname)
---error()
			end

			caidx=groups[ig].caidx
				--- if oxgen atom is set to the residue data
			oa=groups[ig].oatom

			if(oa) then
				span.mid=XYZgetOA(atomxyz,span.smid)  --- try this 
				if(not span.mid) then
					print("span mid error ",span.smid)
				end
				span.midvec,log=XYZvecgetOA(allvec,span.smid, span.mid)   --- default span.mid
			end
			  --- next residue exist
			if(atomxyz[span[2]]) then
				span.ts =XYZcalcavg(atomxyz,span.isp1,span.isp2)
				span.tsvec,log =XYZveccalcavg(allvec,span.isp1,span.isp2, nil) ---default span.ts)
				span.te =XYZcalcavg(atomxyz,span.iep1,span.iep2)
				span.tevec,log2 =XYZveccalcavg(allvec,span.iep1,span.iep2, span.te)  --- default
			else
								print("ISP2",span.isp2,span.isp1)
				span.ts =XYZcalcavg(atomxyz,span.isp1,span.isp1)
				span.tsvec,log =XYZveccalcavg(allvec,span.isp1,span.isp1, nil ) ---default span.ts)
				c_atom=groups[ig].c_atom
				if(c_atom) then
					span.te=XYZgetC(atomxyz,span.smid) 
					--- TEST
					if(not span.te[span.smid]) then print("C SEQ err",span.smid, atomxyz[span.smid]) end

					span.tevec,log=XYZvecgetC(allvec,span.smid,  span.te) -- default
				else
					print("error no c atom ",ig,span.smid)
				end

			end
			if(DEBUG) then
				print("resxyz",span.mid[1], span.mid[2], span.mid[3])
				print(span.ts,span.te)
			end
		  --- for multi residue group
		else
			span.gvec,log =XYZveccalcavg(allvec,span[1],span[2])
			if(not span.grav) then
				print("warning this pdb file has no atoms for group ",ig,span[1],span[2],span.cname)
					--- to stop the job
				error()
			end

			span.gvec,log =XYZveccalcavg(allvec,span[1],span[2])

			 -- print("D",groups[ig].caidx)
			caidx=groups[ig].caidx
				--- if oxgen atom is set to the residue data
			oa=groups[ig].oatom

			if(oa) then
				span.mid=XYZgetOA(atomxyz,span.smid)  --- try this
				span.midvec,log=XYZvecgetOA(allvec,span.smid, span.mid)
			else
				  --- ts,te,mid
				print("no oatom",span[1])
				span.mid=XYZcalcavg(atomxyz,span.smid,span.smid)
				span.midvec,log =XYZveccalcavg(allvec,span.smid,span.smid)
			end


			span.ts =XYZcalcavg(atomxyz,span.isp1,span.isp2)
			span.tsvec,log =XYZveccalcavg(allvec,span.isp1,span.isp2)


			span.te =XYZcalcavg(atomxyz,span.iep1,span.iep2)

			if(not span.te) then
				print("NO te",span.iep1,span.iep2,"mid",span.smid)
			end
			span.tevec,log2 =XYZveccalcavg(allvec,span.iep1,span.iep2)
			 -- print("tsvec",span.ts, span.tsvec)
		end

							if(DEBUG) then 
			 			print("tsvec",span.ts[1], span.tsvec[1][1])
			 			print("tevec",span.te[1], span.tevec[1][1])
						print("LOG",log, "LOG2",log2) print(span.sst,span.isp1,span.isp2,  span.iep1,span.iep2, span.smid,oa)
						print("ts",allvec[1],  span.tevec[1][1])
						print("mid",span.mid[1],span.midvec[1][1])
					end
	end

	--- remove empty grous
	log=nil
	for ig=getn(groups),1,-1 do
		span=groups[ig]
		if(not span.grav) then
			if(not log) then 
				log=""
			end
			log=log..string.format(" %d",ig)
			table.remove(groups,ig)
		end
	end
	if(log) then
		print("delete empty groups in DSSP "..log)
	end

	  --- make bone data to be saved as a BVH 
	bonedata=makeBonedataGroups(groups)

	print("het",#ligand)
	boneadd=mkBoneAddHet(ligand)
	--------------
		--- set the motion data for the bone
	  --- INPUT  :the 3 atom trajectory lists for pt1,pt2,pt3
	  --- OUTPUT :bonedata[].motions[]
	makeBoneTraj(bonedata)
	bonedata.name=PDBNAME
	writeBoneTrajBVH(bonedata,boneadd)
end

  --- support multi chain pdb data
function MakeResgroups(mol)
				---	group is an array of list: {[1]=from, [2]=to, sst=, name=,cname=,term=}
	local resgroup={}
	local groupone
			  -- residue number index resxyz[resno] --- a sparse array
	ir=0
	ic=""
	for ia=1,#mol do
		atom=mol[ia]
		if(not atom.het) then
			resno=mol[ia].resno
			cname=mol[ia].cname
			 -- print(ia,resno,mol[ia].aname,resno,ir)
			if(resno~=ir) or (cname~=ic) then
					--- check if contig
				if(resno==ir+1) and (cname==ic) then
					groupone.pep=resno
					  --- write into the previous residue table
				end
				ir=resno
				ic=cname
				cname=mol[ia].cname or "@"
				  --- set contig flag to the previous group
				groupone={ir,ir+1,sst="p",name=cname..ir,cname=cname}
				table.insert(resgroup,groupone)
									if(ir>150) and (ir<170) then print(#resgroup,groupone,ir,cname,mol[ia].resname) end

			end

			if(mol[ia].aname==" CA ") then
				groupone.caidx=ia

			elseif(mol[ia].aname==" O  ") then
				groupone.oaidx=ia  --- index of the first pdb data
				groupone.oatom=mol[ia]
			elseif(mol[ia].aname==" C  ") then
				groupone.cidx=ia  --- index of the first pdb data
				groupone.c_atom=mol[ia]

			end
		end
	end

	--- debug
	print("residue groups ",ir)
	return resgroup
end

--------------------

  --- rewrite residue number that is the same as reference pdb
  --- comparering residue name
  --- if the residue is not found, take them as biological unit replicated
  --- and  put the new residue number

function AlignPdbSeqXXX(refpdb,mol)

	unitseg=0
	log=""
	ia=1
	oldresno=0
	print("debug",getn(refpdb),#refpdb)
	print(refpdb[1].resname)
	natom=getn(mol)

	rcount=0
	isd=1
	resno=0
	for ia=1,natom do
		xresno=mol[ia].resno
		xresname=mol[ia].resname

		if(oldresno==xresno) then	--- repeat the same residue
			mol[ia].resno=resno
			mol[ia].oldresno=oldresno
		else
			--- find correct residue number in referene pdb
			--- starting from previous matched residue number
			while(true) do
				if(resno==refpdb[isd].resno) then  --- skip this residue
					isd=isd+1
				else
					resname=refpdb[isd].resname
					if(resname==xresname) then
						resno=refpdb[isd].resno
						rcount=rcount+1
						 --- print("debug",ia,isd,resno,resname)
						break
					else
						isd=isd+1
					end
				end
				if(not refpdb[isd]) then
					print("refpdb error final resno", refpdb[isd-1].resno,mol[ia].resname,mol[ia].resno)
					  --- we may stop here,
					  --- but we can still try with duplicated pdb as multimer
					isd=1
					  --- assign different residue number for this unit
					unitseg=unitseg+math.ceil(resno/1000)*1000
					nac1=mol[ia].resno
					res1=mol[ia].resname
					print("residue not found. try repeating the same unit",unitseg)
					print(log) log=""
				end
			end
			 --- set new residue number
			oldresno=mol[ia].resno
			mol[ia].resno=resno+unitseg
			mol[ia].oldresno=oldresno
			log=log..string.format("%d%s ",resno+unitseg,resname )
		end
	end

	print(log)
	print("residue number check ",oldresno,resno,rcount)

	  --- if another chain is added to PDB as a biological assembly, maybe dimmer

	if(unitseg>0) then

			print("try updated numbering for the dimer model",nac1,res1,natom)
			newchain="X"

			print("new residue numbers are asigned starting from  ",unitseg)
			BAssembly={
				newchain="X",
				resnoadd=unitseg,

			}
			print(" XXX numbers of atom may not match reference PDB data.",nac1,natom)
			--return false
		end
	--- ok success
	return true
end

function mkBoneAddHet(ligand)
	bonesadd={}
	newhet=nil
	newchain=nil
	for ia=1,#ligand do

		atom=ligand[ia]
		if ligand[ia].resno==newhet and atom.cname==newchain then
			newbone.to=ia
			  -- use an atom element C 
			if( string.sub(atom.aname,2,2)=="C") then
				newbone.bs={x=atom.x-newbone.pt1.x, y=atom.y-newbone.pt1.y, z=atom.z-newbone.pt1.z}
			end
		else
			newhet=atom.resno
			newchain=atom.cname
			newbone={from=ia,
				to=ia,
				name="h-"..atom.cname..string.format("%d",atom.resno),
				resno=atom.resno,
				pt1={x=atom.x, y=atom.y, z=atom.z},
				bs={x=0,y=0,z=1}
				}
			table.insert(bonesadd,newbone)
			print("hetero atom residue",newhet,newbone.name,atom.cname)
		end
	end
	print("added",#bonesadd)
	return bonesadd
end

------------
  --- calc average xyz atomic coordinate for the residue renge
function XYZgetOA(xyz,pos1)
	if(not xyz[pos1]) then
		  --- print("residue error",pos1)
		return
	end
	oatom=xyz[pos1].oatom
	if(not oatom) then
		--- print("oatom xyz error",pos1)
	end
	return {oatom.x,oatom.y,oatom.z}
end
function XYZvecgetOA(allvec,p1, default)
		local log
		local grav
		  --- for all mode
		log=""
		gvec={}
		for im=1,getn(allvec) do
			grav=XYZgetOA(allvec[im],p1)
			if(not grav) then
				 --- print("no OA")
				grav=default
			end
			table.insert(gvec,grav)
		end
		if(not grav) then
			return nil
		end
		log=log..string.format("%7.3f %7.3f %7.3f ",grav[1],grav[2],grav[3])
		return gvec,log
end

  --- calc average xyz atomic coordinate for the residue renge

  -----------copy and paste hack
function XYZgetC(xyz,pos1)
	local c_atom
	if(not xyz[pos1]) then
		 --- print("C sequence error",pos1,xyz)
		---for idx,cnt in pairs(xyz) do print(idx) end
	else
		c_atom=xyz[pos1].c_atom
	end
	if(not c_atom) then
		 --- print("C atom xyz error",pos1)
		return nil
	end
	return {c_atom.x,c_atom.y,c_atom.z}
end
function XYZvecgetC(allvec,p1, default)
		local log
		local grav
		  --- for all mode
		log=""
		gvec={}
		for im=1,getn(allvec) do
			grav=XYZgetC(allvec[im],p1)
			if(not grav) then
				 --- print("no C")
				grav=default
			end
			table.insert(gvec,grav)
		end
		if(not grav) then
			return nil
		end
		log=log..string.format("%7.3f %7.3f %7.3f ",grav[1],grav[2],grav[3])
		return gvec,log
end

  --- calc average xyz atomic coordinate for the residue renge
function XYZcalcavg(xyz,pos1,pos2)
	local gx,gy,gz = 0,0,0
	local n=0
	  --- 
	for ir=pos1,pos2 do
		 --- 
		if(not xyz[ir]) then
			print("out of range XYZcalcavg()",ir,pos1,pos2)
		elseif(not xyz[ir][3]) then
			print("xyz[3] error",ir,xyz[ir][1],xyz[ir][2],pos1,pos2)
			---error(-1)
		else
				 --- print("DEBUGatm",ir,xyz[ir][1],xyz[ir][2])
			n=n+1
			gx=gx+xyz[ir][1]
			gy=gy+xyz[ir][2]
			gz=gz+xyz[ir][3]
		end
	end
	if(n>0) then
		gx=gx/n
		gy=gy/n
		gz=gz/n
		return {gx,gy,gz}
	end
end

function XYZveccalcavg(allvec,p1,p2,gravdefault)
		local log
		local grav
		  --- for all mode
		log=""
		gvec={}
		for im=1,getn(allvec) do
			if(not allvec[im][p1]) then
				if(im==1) then	print("missing atom",p1,gravdefault) end
				grav=gravdefault
			else
				grav=XYZcalcavg(allvec[im],p1,p2)
								--- print("veccalcavg",p1,p2,grav[1],grav[2],grav[3])
			end
			table.insert(gvec,grav)
		end
		if(not grav) then
			return nil
		end
		log=log..string.format("%7.3f %7.3f %7.3f %d-%d",grav[1],grav[2],grav[3],p1,p2)
		return gvec,log
end
	------ XYZ data util
function XYZadd(a,b,factor)
	if(factor) then
		a.x=a.x+b.x*factor
		a.y=a.y+b.y*factor
		a.z=a.z+b.z*factor
	else
		a.x=a.x+b.x
		a.y=a.y+b.y
		a.z=a.z+b.z
	end
end


-----------------

  --- BHV file write
  --- use  the primary bone: 
  --- head pos at (OFFSET) tail pos at the mid of joints
  ---
bvhf={}

bvhf.h_fmt=[[HIERARCHY
ROOT bvh%s
{
	OFFSET %.3f %.3f %.3f
	CHANNELS 0
]]


--- take ZXY euler in blender (Z first)
---		CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation


  --- joints position relative to the primary node
  --- joint tailer relative to OFFSET
  --- the center of the rotation is on the first OFFSET

  -- try with ZXY
bvhf.b_zxy_fmt=[[	JOINT %s
	{
		OFFSET %.3f %.3f %.3f
		CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
		End Site
		{
			OFFSET %.3f %.3f %.3f
		}
	}
]]
  -- static bone
bvhf.bfix_fmt=[[	JOINT %s
	{
		OFFSET %.3f %.3f %.3f
		End Site
		{
			OFFSET %.3f %.3f %.3f
		}
	}
]]
  --- for connected child bones
bvhf.c_fmt=[[	JOINT %s
	{
		OFFSET %.3f %.3f %.3f
		CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
]]

    --- the trajectory coordinates for nodes of each bones

bvhf.m_fmt=[[MOTION
Frames: %d
Frame Time: %.5f
]]
 ---

function makeBonedataGroups(groups)

	bonedata={
	}

			------------------------------
	  --- make bones from group
	imode=1
	for ig=1,getn(groups) do
		span=groups[ig]

		--print("debug",span[1],span[2])
		--- read the control point mid
		tab={0,0,0,span.mid[1],span.mid[2],span.mid[3]}
		midlist=tab
		midlist[1]=ig
		midlist.xyz=span.mid
		midlist.vec=span.midvec

			--- the start pos
		tab={0,0,0,span.ts[1],span.ts[2],span.ts[3]}
		firstlist=tab
		firstlist.vec=span.tsvec
		firstlist.xyz=span.ts

				---DEBUG
				if(not span.te) then
					print("error NO te", span[1],span[2])
				end
			--- the last pos
		tab={0,0,0,span.te[1],span.te[2],span.te[3]}
		lastlist=tab
		lastlist.vec=span.tevec
		lastlist.xyz=span.te

		name=(span.sst or "L").."-"..span.name
		
		x=midlist.xyz[1]
		y=midlist.xyz[2]
		z=midlist.xyz[3]
			---print("debug",span[1],span[2],span.name)
			--- test
		mx= (lastlist.xyz[1]+firstlist.xyz[1])/2
		my= (lastlist.xyz[2]+firstlist.xyz[2])/2
		mz= (lastlist.xyz[3]+firstlist.xyz[3])/2
		
			--- for debug purpose
			--- print("F",lastlist[4],firstlist[4],x,y,z,mx,my,mz)

		dx= lastlist.xyz[1]-firstlist.xyz[1]
		dy= lastlist.xyz[2]-firstlist.xyz[2]
		dz= lastlist.xyz[3]-firstlist.xyz[3]
		--	boneone={pos={x=tab[1],y=tab[2],z=tab[3]} }
		boneone={
				axis={x=dx,y=dy,z=dz},
				---ctr={x=x,y=y,z=z},
				name=name,
				resno=span[1],
				pep=span.pep  --- peptide link infor
		}

		--- the starting pos
		firstp={x=firstlist.xyz[1],y=firstlist.xyz[2],z= firstlist.xyz[3],debug=span.ts[1]}
		---midp	={x=mx,	y=	my,z=	 mz}
		lastp={x=lastlist.xyz[1],y=lastlist.xyz[2],z= lastlist.xyz[3],debug=span.te[2]}
		
		boneone.span=span --- from/to data
		boneone.pt1=firstp
		boneone.pt2=lastp
		boneone.pt3={x=x,y=y,z=z,debug=0} --- midlist.xyz


---DEBUGif(ig==1) then print("M",x,y,z) print("pt1",firstp.x,firstp.y) end

		---test
		if(span.oatom) then
			  --- print("span oatom",oatom.x,ig)
--- DEBUG boneone.pt3={x=oatom.x, y=oatom.y, z=oatom.z}
			  --- print(firstp.x,firstp.y)
			  --- print(lastp.x,lastp.y)
			  --- print("traj")
			  --- print(firstlist.vec[1][1],firstlist.vec[1][2])
			---print(firstlist.vec[2][1],firstlist.vec[2][2])
			---print(firstlist.vec[3][1],firstlist.vec[3][2])
		end
		-- print(ig,"debug",firstp.debug,lastp.debug)

		--boneone.pt3=span.oatom 	--- JSUT TEST
		  --- if(ig==1) then print("XYZ",x,y,z) end --- DEBUG
		
---TRAJVESION
	---print("DEBUG",ig,span.tsvec)
if(not firstlist.vec) then
	--print("error firstlist.vec",firstlist,boneone.name)
	--print(span.tsvec,span.ts)
end
		firstp.traj=firstlist.vec
		lastp.traj=lastlist.vec
		boneone.pt3.traj=midlist.vec
----------------
		-------------------
		table.insert(bonedata,boneone)

	end

	return bonedata
end

  --- calc tans+rot for each bones of all trajectory
function makeBoneTraj(bonedata)
	local bonetraj
	--- from bone motions as 3point pt1.traj, pt2.traj, pt3.traj in 
	--- in [123] fromat is converted
	--- bonetraj[frame][bone].pt1  .pt2  .pt3 {xyz}format
	--- each bone trajectores are stored at bone.motions[]
	--------------

	ntraj=getn(bonedata[1].pt3.traj)

	for ib=1,#bonedata do
		bonedata[ib].motions={}
	end
	
	local am,mb,sb,mb,eb
	
	local pmop=Mop:new()
	  --- the 1st bonetraj is the initial model

	for it=1,ntraj do
					---t=(it-1)*deltat
		for ib=1,#bonedata do
			  --- calc 3 division displacement vectors
			  --- make composit vectors for the bone
			midp  =bonedata[ib].pt3
			firstp=bonedata[ib].pt1
			lastp =bonedata[ib].pt2

			  --- take trajectory

			if(bonedata[ib].pt1.traj) then
				firstp2=XYZval(bonedata[ib].pt1.traj[it])
			else
				firstp2=nil
			end
			if(bonedata[ib].pt2.traj) then
				lastp2=XYZval(bonedata[ib].pt2.traj[it])
			else
				lastp2=nil
			end
			if(bonedata[ib].pt3.traj) then
				midp2=XYZval(bonedata[ib].pt3.traj[it])
			else
				midp2=nil
			end

			--- TEST coord
			if(DEBUG) then
				print(" DEBUG p2",ib,it,bonedata[ib].name)
				print(firstp2,midp2,lastp2) 
				print(firstp2.debug,midp2.debug,lastp2.debug)
				debug3pt(firstp,midp,lastp,firstp2,midp2,lastp2)
			end

			  --- list bones
			if(DEBUG) then
				len=bonedata[ib].span[2]-bonedata[ib].span[1]+1
				  -- print(ib, len,bonedata[ib].name,bonedata[ib].span[1])
			end

			  --- 
			if(not firstp2) then  --- lack trajectory take the same one as the previous residue
				--rx,ry,rz=0,0,0
				--tx,ty,tz=0,0,0

				tx,ty,tz=mop:calc(firstp.x,firstp.y,firstp.z)

				trajone={x=tx,y=ty,z=tz
					,dx=tx-firstp.x,dy=firstp.y,dz=firstp.z
					,xrot=rx,yrot=ry,zrot=rz}
				if(interpmode) then
					interpmode2=ib
					--- if(it==2) then print(ib,interpmode,bonedata[ib].resno) end
				else
					  --- save the bone motion should be recalculated later
					interpmode=ib
					interpmode2=ib
				end
			else
				if(interpmode) then
					--- if previous bone lacks motion data, try interpolating them
					--- save them in mixing table
					if(it==2) then
						print("residue no motions",bonedata[interpmode].resno,bonedata[interpmode2].resno,
									bonedata[interpmode2].pep)
						  --- pep : peptide bond exist to the next bone
						bonedata[ib].mixing={bstart=interpmode,bend=interpmode2,mop=mop}
					end
					interpmode=nil
					  --- no interpolation needed from this bone
				end
				  --- matrix to fit 3 points
				mop,elog=fit3points(firstp,lastp,midp,firstp2,lastp2,midp2)	 

				if(DEBUG) then 
					  -- print(vlist[im][1],im,bonedata.nmodefreq[im])
					  ---elog=""
					  
					end

				--- if ib<10 then elog=elog or "DEBUG" end
				if(elog) then
					print("fit3p",bonedata[ib].name,it,ib,bonedata[ib].span[1],
					  bonedata[ib].span[2],elog)
					print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
					  firstp.x,firstp.y,firstp.z,midp.x,midp.y,midp.z))
					print(string.format("%8.3f %8.3f %8.3f ",lastp.x,lastp.y,lastp.z))
					print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
					  firstp2.x,firstp2.y,firstp2.z,midp2.x,midp2.y,midp2.z))
					print(string.format("%8.3f %8.3f %8.3f ",lastp2.x,lastp2.y,lastp2.z))
	---error()
				end 

			  	 --- take the location relative to the center of this object
				
				tx=mop.trn[1]
				ty=mop.trn[2]
				tz=mop.trn[3]

				if(FMTB) then
									-- the center of rotation is pt1

					tx=firstp2.x-firstp.x
					ty=firstp2.y-firstp.y
					tz=firstp2.z-firstp.z
				 ---- for a child bone ---
				else
						--- not debuged yet	
						--- compensate rotation of previous bone	
					mop.mtx=YotM3xM3(pmop.mtx,mop.mtx)
					pmop.mtx=YotM3Trans(mop.mtx)
					  --- just debug
					---mop.mtx=EulerMtx(0,0,0)
				end

				   --- note : BVH format takes ZXY rotation in Blender
				   --- first Z then X , last Y
				rx,ry,rz=YotM3GetYxzDeg(mop.mtx)
			 
				if(elog) then
					print("deg",rx,ry,rz)
				end
			
---			  --- for FMTB use global coordinate
				  --- for FMTA use local coordinates of the begining bone
			--tx=tx+bonedata[ib].pt1.x	-- -boneCx
			--ty=ty+bonedata[ib].pt1.y	-- -boneCy
			--tz=tz+bonedata[ib].pt1.z	-- -boneCz

				ex,ey,ez=mop:calc(lastp.x,lastp.y,lastp.z)

			  --- save the result
				trajone={x=tx+bonedata[ib].pt1.x,y=ty+bonedata[ib].pt1.y,z=tz+bonedata[ib].pt1.z
					,dx=tx,dy=ty,dz=tz
					,ex=ex,ey=ey,ez=ez
					,xrot=rx,yrot=ry,zrot=rz}
			end
			bonedata[ib].motions[it]=trajone
			--- the calc mixings for previous bones
			if(bonedata[ib].mixing) then
				mixing=bonedata[ib].mixing
									  --- the bone to follow
				bnext=bonedata[ib]

				if(it>=2) and bnext.motions then
							--- print("mixcalc",mixing.bstart,mixing.bend)
					itb1=bonedata[mixing.bstart]
					itb2=bonedata[mixing.bend]

					bprev=bonedata[mixing.bstart-1]


						--- print("ib",bprev.resno,bnext.resno)
				  --- matrix to fit 3 points
				  firstp =bprev.pt2
				  firstp2=XYZval(bprev.pt2.traj[it])
					lastp  =bnext.pt1
					lastp2 =XYZval(bnext.pt1.traj[it])
					midp   =bnext.pt2
					midp2  =XYZval(bnext.pt2.traj[it])

							--- print("prev",bprev.resno,firstp2.x,firstp2.y,firstp2.z)
					mop,elog=fit3points(firstp,lastp,midp,firstp2,lastp2,midp2)	

					if(it==2) then
						print("interpolate range",itb1.resno,itb2.resno,"tward",bnext.resno)
					end

					rx,ry,rz=YotM3GetYxzDeg(mop.mtx)
					  --- the gap size

					alebone=0
					rgap=bnext.resno-itb1.resno

					tx2=firstp2.x
					ty2=firstp2.y
					tz2=firstp2.z

					tx2=bprev.motions[it].ex
					ty2=bprev.motions[it].ey
					tz2=bprev.motions[it].ez

					for jb=mixing.bstart,mixing.bend do
						bonedata[jb].motions[it].x=tx2 --- +bonedata[jb].pt1.x
						bonedata[jb].motions[it].y=ty2 --- +bonedata[jb].pt1.y
						bonedata[jb].motions[it].z=tz2 --- +bonedata[jb].pt1.z
								--- print("tx",tx2,ty2,tz2)
					 	bonedata[jb].motions[it].xrot=rx
					 	bonedata[jb].motions[it].yrot=ry
					 	bonedata[jb].motions[it].zrot=rz
						  --- mixing alpha for start bone and end bone
						alsbone=alebone
						alebone=(bonedata[jb].resno-itb1.resno+1)/rgap

							-- the start pos firstp2 
							-- calc the end pos to use following bone
						secondp=bonedata[jb].pt2
						tx2,ty2,tz2=mop:calc(secondp.x,secondp.y,secondp.z)



									if(false) then
										tx2=tx2- secondp.x+firstp.x
										ty2=ty2- secondp.y+firstp.y
										tz2=tz2- secondp.z+firstp.z

										bnext=bonedata[jb+1]  --- the next bone 
										nbmot=bnext.motions[it]
										tx2=nbmot.x- secondp.x+firstp.x
										ty2=nbmot.y- secondp.y+firstp.y
										tz2=nbmot.z- secondp.z+firstp.z
									end
						--- print("motion estimated",ib,jb,"res",bonedata[jb].resno) 

					end

				end
			end
---
		end
	end
	return 
end

function debug3pt(firstp,midp,lastp,firstp2,midp2,lastp2,mop)
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	  midp.x,midp.y,midp.z,midp2.x,midp2.y,midp2.z))
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	  firstp.x,firstp.y,firstp.z,firstp2.x,firstp2.y,firstp2.z))
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	  lastp.x,lastp.y,lastp.z,lastp2.x,lastp2.y,lastp2.z))
	if(mop) then
		print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
		mop.mtx[1], mop.mtx[2], mop.mtx[3], mop.mtx[4], mop.mtx[5], mop.mtx[6]))
		print(string.format("%8.3f %8.3f %8.3f",
		mop.mtx[7], mop.mtx[8], mop.mtx[9] ))
	end
end


--- On importing bvh in blender one root bone will be added at
--- the center of all bones.
--- Adding the size of the root bone here will locate all bones
--- correct location

BONELOCx=2
BONELOCy=0
BONELOCz=0

function PatchBoneLoc(x,y,z)
  x=x+BONELOCx
  y=y+BONELOCy
  z=z+BONELOCz
  return x,y,z
end

FMTB=1
  --- save the bone motions data to a BVH file
--------------
function writeBoneTrajBVH(bonedata,boneadd)
	--- local boneCx,boneCy,boneCz,bonec
		--- calc center of all bone
	boneadd=boneadd or {}

	nbones=#bonedata
	bonec={x=0,y=0,z=0}
	for ib=1,nbones do
		XYZadd(bonec,bonedata[ib].pt3)
	end

	boneCx=bonec.x/nbones
	boneCy=bonec.y/nbones
	boneCz=bonec.z/nbones
	  --- patch the bone coordinate offset from the
	  --- center of all bone 

	boneCx,boneCy,boneCz=PatchBoneLoc(boneCx,boneCy,boneCz)

	bonedata.pos={x=boneCx,y=boneCy,z=boneCz}

	print("center pos for the root bone",boneCx,boneCy,boneCz)

	  --- just incase the vector is nil
	noaxis={x=5,y=0,z=0}

		-----------
		
	nframes=getn(bonedata[1].motions)
	tframe=0.04 --- a tentative value, the time step 1/25
	nbones=#bonedata
	io.write(string.format(bvhf.h_fmt,PDBNAME,
			 boneCx,boneCy,boneCz))		--- NO SWAP ZY

	
	ox=0
	oy=0
	oz=0
	bclose=0
	for ib=1,nbones do
		name=bonedata[ib].name
		b1=bonedata[ib].pt1 --- pt3			---TEST
		
		bs=bonedata[ib].axis or noaxis

	  FMTB=1  --- force all bone independent 
		if(FMTB) then
			io.write(string.format(bvhf.b_zxy_fmt,name,
				b1.x-boneCx,b1.y-boneCy,(b1.z-boneCz),
				bs.x,bs.y,bs.z))						
		else
			x=b1.x-boneCx
			y=b1.y-boneCy
			z=b1.z-boneCz
			io.write(string.format(bvhf.c_fmt,name,x-ox,y-oy,z-oz))
			ox=x
			oy=y
			oz=z
	
			bclose=bclose+1

		end				
	end
	if(bclose>0) then
		  --- the last bone size
		ox=bs.x
		oy=bs.y
		oz=bs.z
		io.write(string.format("		End Site\n    {\n    	OFFSET %.3f %.3f %.3f\n    }\n",ox,oy,oz))

		io.write(string.rep("}",bclose))
		io.write("\n")

	end
	  --- additional static bones
	for ib=1,#boneadd do
		name=boneadd[ib].name
		b1=boneadd[ib].pt1 
		bs=boneadd[ib].bs or noaxis

		io.write(string.format(bvhf.bfix_fmt,name,
				b1.x-boneCx,b1.y-boneCy,(b1.z-boneCz),
				bs.x,bs.y,bs.z))
	end


	io.write("}\n") --- end of skelton
	
	io.write(string.format(bvhf.m_fmt,nframes+1 ,tframe))
	for im= 0 ,nframes do  --- frame 0 initial model is added 

		ox=boneCx
		oy=boneCy
		oz=boneCz
		for ib=1,nbones do
			locrot=bonedata[ib].motions[im]
			if(not locrot) then
				--- print("ERR",im,ib)  --- for the first frame with the rest pose 
				locrot=bonedata[ib].pt1
				--locrot={x=bonedata[ib].x,y=0,z=0,xrot=0,yrot=0,zrot=0}
			end

			if(FMTB) then
				io.write(string.format("%.5f %.5f %.5f",locrot.x-boneCx,locrot.y-boneCy,locrot.z-boneCz))
			else

				  --- relative coordinate
				io.write(string.format("%.5f %.5f %.5f",locrot.x-ox,locrot.y-oy,locrot.z-oz))
				ox=locrot.x
				oy=locrot.y
				oz=locrot.z

			end

			xrot=locrot.xrot or 0
			yrot=locrot.yrot or 0
			zrot=locrot.zrot or 0
-----
			io.write(string.format(" %.3f %.3f %.3f ",zrot,xrot,yrot))

		end
		io.write("\n") --- end of skelton
	end
	io.write("\n") --- end of skelton

end

function XYZval(ary)
	if(ary) then
		return {x=ary[1],y=ary[2],z=ary[3]}
	end
end

--------------------
main()
 
