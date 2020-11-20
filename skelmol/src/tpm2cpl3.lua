 
  --- make a trajectory  motion bones from pdb trajectory 
  ---
  --- the updated coord n Lua by 2018 April 11 Y.Ueno
  --- the coord originated by 2014 Oct24 Yutaka Ueno and Shinya Muraoka
  --- AIST and NAIST Japan

  --- a sample code
  --- bvh-bone data generator 
  --- this version takes DSSP characterization for residues 
  --- to make groups of peptide in the protein
  --- ProModeElastic results are on the web http://pdbj.org/promode-elastic
  --- provided by Wako&Endo, Waseda U.


    
dofile("fit3p.lua")
dofile("pdbmol.lua")

  --- some compat stuff
function getn(t) if(t) then return #t end end


function main()
	print(arg,argv)
	print(arg[0],arg[1])
---if(1) then return end

	molname="1amy"

	--- first get DSSP
	DFILE=molname.."-dssp.txt"
	--- but this time take another one
	  --- get 2nd st DSSP notation 

	TFILE=""..molname.."0_md300.pdb"
	BFILE="tpm-"..molname..".bvh"
	BFILE="tpm-"..molname.."0_md300.bvh"

	io.output(BFILE)

	io.input(TFILE)

			 --- only CA atoms
	allpdb={}
	repeat
		mol=msReadPdb()
		--- get only CA
		for ia=getn(mol),1,-1 do
			aname=mol[ia].aname
			if(aname==" CA ") then
				-- ok
			else
				table.remove(mol,ia)
			end
		end

		if(mol.modelno) then
			table.insert(allpdb,mol)
			--- load next MODEL data
			print("MODEL",mol.modelno,#mol)
		else
			mol=nil
			break
		end
	until not mol

------
	  --- numbers of residues

--- convert xyz to array
	--- calc movements vector for each traj
	nframe=getn(allpdb)
	mol=allpdb[1]
	nra=getn(mol)
print("NRES",nra,"models",nframe)
	allvec={}
	for is=1,nframe do
		vec={}
		table.insert(allvec,vec)
		mol2=allpdb[is]
		for ia=1,nra do
			atomone=mol[ia]
			atomtwo=mol2[ia]
			xyz={}
			table.insert(vec,xyz)
			if(is==1) then
				xyz[1]=atomone.x
				xyz[2]=atomone.y
				xyz[3]=atomone.z
			else
				xyz[1]=atomtwo.x -- -atomone.x
				xyz[2]=atomtwo.y -- -atomone.y
				xyz[3]=atomtwo.z -- -atomone.z
			end
		end
	end

	  --- set xyz trajectory 
	atomxyz=allvec[1]
	table.remove(allvec,1)

--------------------
	  --- 2nd st DSSP notations are in the PROMODE
	  --- So, simply provide CA atom list for the residue list
	dssplist=ParseDSSP(DFILE)
	segres=GetDSSPspan(dssplist) --- make a segment list 
	  --- create group list from the DSSP
	groups=MakeDSSPgroups(segres)


--------------------same as ela2cpl3

	  --- calc gravity and vector for groups from sst 
	bonelist={}
	for ig=1,getn(groups) do
		span=groups[ig]
		  --- set 1st, last, and mid-control points
		  --- for helix or sheet, apply long range bone
		  --- the control point smid
		smid=math.floor((span[2]+span[1])/2)

		-- special for the alpha helix, an emprical tip
		---if(sst=="H")or(sst=="E") then
		if(span.sst=="H") then
			smid=span[1]+5
		else
			smid=span[1]+1 ---4
		end
----trial
----		if(span.sst=="S") then
----			smid=span[1]+1
----		end
		---
		len=span[2]-span[1]+1
		isp1=span[1]
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
		if(span.sst~="H")and(span.sst~="E") then
			isp2=smid-1
			iep1=smid+1
		end

		span.smid=smid
		span.isp1=isp1
		span.isp2=isp2
		span.iep1=iep1
		span.iep2=iep2
	end

	for ig=1,getn(groups) do
		span=groups[ig]

		  --- set bone coordinate from atom xyz
		span.grav=XYZcalcavg(atomxyz,span[1],span[2])
		span.gvec,log =XYZveccalcavg(allvec,span[1],span[2])

		  --- ts,te,mid
		span.mid=XYZcalcavg(atomxyz,span.smid,span.smid)
		span.midvec,log =XYZveccalcavg(allvec,span.smid,span.smid)

		span.ts =XYZcalcavg(atomxyz,span.isp1,span.isp2)
		span.tsvec =XYZveccalcavg(allvec,span.isp1,span.isp2)

		span.te =XYZcalcavg(atomxyz,span.iep1,span.iep2)
		span.tevec =XYZveccalcavg(allvec,span.iep1,span.iep2)
	end

	group2bvh(groups,freqlist,molname)
end


  --- make a span list from residue list, alpha carbon list
function GetDSSPspan(reslist)
	
	glist={}
	inthem=""
	gno=0
	nra=getn(reslist)
	for ia=1,nra do
		sst=reslist[ia].sst
		rno=reslist[ia].rno
		if(sst~=inthem) then
			inthem=sst
			gno=gno+1
			cname=reslist[ia].cname
			glist[gno]={num=0,sst=sst,from=ia,to=ia, cname=cname }
		end
		---reslist[ia].gno=gno
		glist[gno].num=glist[gno].num+1
		glist[gno].to=ia
	end
	inchain=""
	for ip=1,getn(glist) do
		if(inchain~=glist[gno].cname) then
			glist[gno].term="n"
			inchain=glist[gno].cname
		end
	end
	return glist
end
function ParseDSSP(dsspfname)
	local dssplist={}
	local fp=io.open(dsspfname,"r")
	local idx,sst,cname,rcode
	if( not fp ) then print("xxx no file",name) return end

	while(true) do
		line=fp:read()    --- 1 coment 
		if(not line) then break end
		key=string.upper(string.sub(line,1,15)) 
		if(key=="  #  RESIDUE AA")then
			break
		end
	end
	while(line) do
		line=fp:read()    --- 1 coment 
		if(not line) then break end
		idx= tonumber(string.sub(line,1,5)) 
		if(idx) then
			rno= tonumber(string.sub(line,6,10)) 
			sst= string.sub(line,17,17)  --- s-struc assignment
			cname= string.sub(line,12,12) -- chain name
			rcode= string.sub(line, 9, 9) -- residue 1char
			-- 3truns, 4turns, 5turns, bend, chirality, betaKL
			-- beta bridge partner, sheet, solvent
			table.insert(dssplist,{rno=rno,sst=sst,cname=cname,rcode=rcode})
		end
		--print(idx,rno,sst)
	end
	fp:close()
	if(sst) then
		return dssplist
	end
end

	--- create groups from DSSP residue list 
function MakeDSSPgroups(segres)
	  --- make DSSP groups of residues
	local groups={}
	sstrec=nil  --- the current group 
	inchain=nil
	for ip=1,getn(segres) do
		ginf=segres[ip]
		sst=ginf.sst
		len=ginf.to-ginf.from
		if(sst=="H")or(sst=="E") then
			sstrec=nil --- new
		  ---turn=="T" bend=="S" 3helix=="G"  5helix=="I"
		else
			if(sstrec) then
				len=ginf.to-sstrec[1]
				  --- just omit special amino acid annotation
				  --- by few residues.. to look into segments 
				if(len>2) then
					sstrec=nil --- new
				end
			end
		end
		if(sstrec) then
			  --- just extend current record
			sstrec[2]=ginf.to
		else
			  --- make a new group
			sstrec={ginf.from, ginf.to,
				sst=sst,
				name=cname .. ginf.from,
				term=ginf.term}
			table.insert(groups,sstrec)
		end
	end

	--- for ig=1,getn(groups) do print(ig,groups[ig].term) end

	  --- this parameter control how many bones with sizes of
	  --- residues for making groups

	--- limit in number of residues
	  --- maximum 
MAXLEN=15
	  --- minimum, these group will be testted for the marge
MINLEN=6
	 --- remove small groups
	 --- except N,C,terminus groups

	 --- first, dissolve small beta-sheet or helix
	 --- and extend loop regions to include them
	for ig=getn(groups),1,-1 do
		span=groups[ig]
		len=span[2]-span[1]-1
		if(len<MINLEN) then
	 		 --- print("SMALL",ig,span.name,len,span.sst,span[2]-span[1])
			  --- marge this into the previous group
			if(ig>1)and(span.term~="n") then
				sst=groups[ig-1].sst
				if(sst=="H")or(sst=="E") then
					--- no marge
				else
					ia=span[2]
					sst=span.sst
					span=groups[ig-1]
					span[2]=ia
					 -- print("marged to",span.name,span.sst)
					if(span.sst==" ") then
						span.sst=sst
					end
				end
				  -- otherwize, we ignore this span
			end
			table.remove(groups,ig)  --- just cut out
		end
	end

	  --- divide very large groups
	for ig=getn(groups),1,-1 do
		span=groups[ig]
		len=span[2]-span[1]

		  ---	if(span.sst=="H")or(span.sst=="E") then
		if(span.sst~="H") then
			if(len>=MINLEN) then
				if(len>MAXLEN) then
					print("seperate",span[3] or "-",span[1],span[2],len,mid)
					  --- divide two
					smid=math.floor((span[2]+span[1])/2)
					span2={smid+1,span[2],"+"}
					span2.sst="+"
					span[2]=smid
					span2.name=span.name
					table.insert(groups,ig+1,span2)
				end
				  --- keep N,C-terminals
				  --- but cut out if small
			elseif(span.term) then
				if(len<2) then
					print("Delete",span.name,span[1],len)
					table.remove(groups,ig)
				end
				
			elseif(not span.term) then
				print("delete",span.name,span[1],span[2]-span[1])
				table.remove(groups,ig)
			end
		end

	end
	return groups
end
-------------

  --- calc average xyz atomic coordinate for the residue renge
function XYZcalcavg(xyz,pos1,pos2)
	local gx,gy,gz = 0,0,0
	local n=pos2-pos1+1
	for ia=pos1,pos2 do
		gx=gx+xyz[ia][1]
		gy=gy+xyz[ia][2]
		gz=gz+xyz[ia][3]
	end
	gx=gx/n
	gy=gy/n
	gz=gz/n
	return {gx,gy,gz}
end

function XYZveccalcavg(allvec,p1,p2)
		local log
		local grav
		  --- for all mode
		log=""
		gvec={}
		for im=1,getn(allvec) do
			grav=XYZcalcavg(allvec[im],p1,p2)
			table.insert(gvec,grav)
		end
		log=log..string.format("%7.3f %7.3f %7.3f ",grav[1],grav[2],grav[3])
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
	OFFSET %.1f %.1f %.1f
	CHANNELS 0
]]

  --- joints position relative to the primary node
  --- joint tailer relative to OFFSET
  --- the center of the rotation is on the first OFFSET
bvhf.b_fmt=[[	JOINT %s
	{
		OFFSET %.3f %.3f %.3f
		CHANNELS 6 Xposition Yposition Zposition Xrotation Yrotation Zrotation
		End Site
		{
			OFFSET %.3f %.3f %.3f
		}
	}
]]

    --- the trajectory coordinates for nodes of each bones

bvhf.m_fmt=[[MOTION
Frames: %d
Frame Time: %.5f
]]
 ---

--- for blender, the size of the parent bone
--- the positoin is automatically calculated by blender
--- as the center of all bones

BONELOCx=2
BONELOCy=0
BONELOCz=0

function PatchBoneLoc(x,y,z)
  x=x+BONELOCx
  y=y+BONELOCy
  z=z+BONELOCz
  return x,y,z
end

function group2bvh(groups,freqlist,molname)

	nmode=getn(freqlist)
	bonedata={
	nmodefreq=freqlist
	}

			------------------------------
	  --- make bones from group
	imode=1
	for ig=1,getn(groups) do
		span=groups[ig]

		
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
				ctr={x=x,y=y,z=z},
		}

		--- the starting pos
		firstp={x=firstlist.xyz[1],y=firstlist.xyz[2],z= firstlist.xyz[3]}
		---midp	={x=mx,	y=	my,z=	 mz}
		lastp={x=lastlist.xyz[1],y=lastlist.xyz[2],z= lastlist.xyz[3]}
		
		boneone.span=span --- from/to data
		boneone.pt1=firstp
		boneone.pt2=lastp
		boneone.pt3={x=x,y=y,z=z} --- midp
---		boneone.pt3={x=mx,y=my,z=mz} --- midp
		
---TRAJVESION

		firstp.traj=firstlist.vec
		lastp.traj=lastlist.vec
		boneone.pt3.traj=midlist.vec
----------------
if(TRAJVERSION) then
		--- store motion data /translation and rotation
		local vlist={}
		local firstv2,midv2,lastv2
		tab=midlist
		
		for ig=1,nmode do
			it=ig*3+4	--- the start
			--	vlist[ig]={tab[it],tab[it+1],tab[it+2]}
--print("debug",ig,firstlist.vec[ig])
			firstv2={x=firstlist.vec[ig][1],y=firstlist.vec[ig][2],z=firstlist.vec[ig][3]}
			midv2	={x=midlist.vec[ig][1],	y=midlist.vec[ig][2],	z=midlist.vec[ig][3]}
			lastv2 ={x=lastlist.vec[ig][1], y=lastlist.vec[ig][2], z=lastlist.vec[ig][3]}

		--- movements for 3 divisions
			vlist[ig]={}
			vlist[ig].vt3=midv2
			vlist[ig].vt1=firstv2
			vlist[ig].vt2=lastv2

			imode=imode+1
		end
end
		-------------------
		boneone.name="sm-"..name
		boneone.nmodevect=vlist
		table.insert(bonedata,boneone)
		
	end
	
function XYZval(ary)
	return {x=ary[1],y=ary[2],z=ary[3]}
end
	--------------
	  --- make motion data 
	  --- INPUT  group[].tsvec esvec midvec 
	  --- OUTPUT:bonedata[].motions[]
	---maketrajNM(bonedata,nframes)
	ntraj=getn(firstlist.vec)
	bonetraj={}
	for it=1,ntraj do
		bonedata2={}
		table.insert(bonetraj,bonedata2)
		for ib=1,getn(bonedata) do
			  --- change array 
			bonetwo={}
			table.insert(bonedata2,bonetwo)
			boneone=bonedata[ib]
			bonetwo.pt1=XYZval(boneone.pt1.traj[it])
			bonetwo.pt2=XYZval(boneone.pt2.traj[it])
			bonetwo.pt3=XYZval(boneone.pt3.traj[it])
		end
	end

	makeBoneTraj(bonedata,bonetraj)

	--------------
		--- calc center of all bone
	nbones=#bonedata
	bonec={x=0,y=0,z=0}
	for ib=1,nbones do
		XYZadd(bonec,bonedata[ib].pt3)
	end

	BoneCx=bonec.x/nbones
	BoneCy=bonec.y/nbones
	BoneCz=bonec.z/nbones
	  --- patch the bone coordinate offset from the
	  --- center of all bone 

	BoneCx,BoneCy,BoneCz=PatchBoneLoc(BoneCx,BoneCy,BoneCz)

	bonedata.pos={x=BoneCx,y=BoneCy,z=BoneCz}

	print("center pos",BoneCx,BoneCy,BoneCz)


	  --- just incase the vector is nil
	noaxis={x=5,y=0,z=0}

		-----------
		
	nframes=getn(bonetraj)
	tframe=0.014 --- a tentative value, the time step
	nbones=#bonedata
	io.write(string.format(bvhf.h_fmt,molname,
			BoneCx,BoneCz,-BoneCy))		--- SWAP ZY
	
	
	--- BoneCx,BoneCy,BoneCz=0,0,0
	
	for ib=1,nbones do
		name=bonedata[ib].name
		b1=bonedata[ib].pt1 --- pt3			---TEST
		
		bs=bonedata[ib].axis or noaxis
		io.write(string.format(bvhf.b_fmt,name,
				b1.x-BoneCx,b1.z-BoneCz,-(b1.y-BoneCy),
				bs.x,bs.z,-bs.y))		 --- SWAP ZY				
				
			--	b1.x-BoneCx-bs.x/2,b1.z-BoneCz-bs.z/2,-(b1.y-BoneCy-bs.y/2),
			--	bs.x,bs.z,-bs.y))		 --- SWAP ZY
	end
	io.write("}\n") --- end of skelton
	
	io.write(string.format(bvhf.m_fmt,nframes,tframe))
	for im=1,nframes-1 do
		for ib=1,nbones do
			locrot=bonedata[ib].motions[im]
			if(not locrot) then
				print("ERR",im,ib)
			end
			io.write(string.format("%.5f %.5f %.5f",locrot.x-BoneCx,locrot.z-BoneCz,-(locrot.y-BoneCy)))
																						--- SWAP ZY
			xrot=locrot.xrot or 0
			yrot=locrot.yrot or 0
			zrot=locrot.zrot or 0

-----
if(ib==159999) then print("FUNY",zrot,im)
print(b1.x,b1.y,b1.z)
print(locrot.x,locrot.y,locrot.z)
end

			io.write(string.format(" %.3f %.3f %.3f ",xrot,zrot,-yrot))
		end
		io.write("\n") --- end of skelton
	end
	io.write("\n") --- end of skelton

end

function makeBoneTraj(bonedata,bonetraj)

	for ib=1,#bonedata do
		bonedata[ib].motions={}
	end
	
	local am,mb,sb,mb,eb
	
	  --- the 1st bonetraj is the initial model
	for it=1,getn(bonetraj)-1 do
		bone2=bonetraj[it+1]
					---t=it*deltat
		for ib=1,#bonedata do
			  --- calc 3 division displacement vectors
			  --- make composit vectors for the bone
			midp  =bonedata[ib].pt3
			firstp=bonedata[ib].pt1
			lastp =bonedata[ib].pt2

			  --- take trajectory
			midp2=bone2[ib].pt3
			firstp2=bone2[ib].pt1
			lastp2=bone2[ib].pt2

			  --- list bones
			if(it==1) then
				len=bonedata[ib].span[2]-bonedata[ib].span[1]+1
				print(ib, len,bonedata[ib].name,bonedata[ib].span[1])
			end
			  --- matrix to fit 3 points
			mop,elog=fit3points(firstp,lastp,midp,firstp2,lastp2,midp2)	 
 
			if(elog) then
				print(bonedata[ib].name,it,ib,bonedata[ib].span[1],
				  bonedata[ib].span[2],elog)
				print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
				  firstp2.x,firstp2.y,firstp2.z,midp2.x,midp2.y,midp2.z))
				print(string.format("%8.3f %8.3f %8.3f ",lastp2.x,lastp2.y,lastp2.z))
			end 
				if(it==1) and (ib==1) then 
				  -- print(vlist[im][1],im,bonedata.nmodefreq[im])
				end
		  	 --- take the location relative to the center of this object
			
			tx=mop.trn[1]
			ty=mop.trn[2]
			tz=mop.trn[3]
			
			   --- note : BVH format takes ZYX rotation
			rx,ry,rz=YotM3GetZyxDeg(mop.mtx)
		 
			if(elog) then
				print(math.deg(rx),math.deg(ry),math.deg(rz))
			end
		
			  --- use local coordinates of the begining bone
			tx=tx+bonedata[ib].pt1.x	-- -BoneCx
			ty=ty+bonedata[ib].pt1.y	-- -BoneCy
			tz=tz+bonedata[ib].pt1.z	-- -BoneCz

			  --- save the result
			trajone={x=tx,y=ty,z=tz,xrot=rx,yrot=ry,zrot=rz}
			bonedata[ib].motions[it]=trajone
---DEBUG
			if(ib==1599) then print("funny",rx,ry,rz,it)
			print(bonedata[ib].name,it,ib,bonedata[ib].span[1],
			  bonedata[ib].span[2],elog)
				debug3pt(firstp,midp,lastp,firstp2,midp2,lastp2,mop)
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
	print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",
	mop.mtx[1], mop.mtx[2], mop.mtx[3], mop.mtx[4], mop.mtx[5], mop.mtx[6]))
	print(string.format("%8.3f %8.3f %8.3f",
	mop.mtx[7], mop.mtx[8], mop.mtx[9] ))
end

--------------------
main()
 
