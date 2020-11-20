 
  --- make a normal mode motion bones from pdb-elastic normal mode 
  ---
  --- the updated coord n Lua by 2017 Oct24 Y.Ueno
  --- the coord originated by 2014 Oct24 Yutaka Ueno and Shinya Muraoka
  --- AIST and NAIST Japan

  --- a sample code
  --- bvh-bone data generator 
  --- this version takes DSSP characterization for residues 
  --- to make groups of peptide in the protein
  --- ProModeElastic results are on the web http://pdbj.org/promode-elastic
  --- provided by Wako&Endo, Waseda U.


  --- the primary node 
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

 ---
    
dofile("splitjoin.lua")
dofile("fit3p.lua")

function main()
print(arg,argv)
print(arg[0],arg[1])
---if(1) then return end

	  --- read 10 files in PDBElastic calc resut tar.gz archive
		---AFILE="1amy-results/1amy_#.flcatmE"
		---AFILE="1a00-results/1a00_#.flcatmE"

	molname="1a00"
	molname="1amy"

	AFILE=molname.."-results/"..molname.."_#.flcatmE"
	BFILE="ela3-"..molname..".bvh"

	atomvec,atomxyz=readNMfile(AFILE ,1,0)
			 --- only CA atoms

	io.output(BFILE)

	  --- numbers of residues
	nra=table.getn(atomxyz)

	freqlist={}
	allvec={}
	table.insert(allvec,atomvec)
	freqlist[1]=atomvec.freq
	for im=2,10 do
		atomvec=readNMfile(AFILE ,im,nra)
		table.insert(allvec,atomvec)
		freqlist[im]=atomvec.freq
	end

	  --- get 2nd st DSSP notation 
	
	glist={}
	inthem=""
	gno=0
	for ia=1,nra do
		sst=atomxyz[ia].sst
		rno=atomxyz[ia].rno
		if(sst~=inthem) then
			inthem=sst
			gno=gno+1
			cname=atomxyz[ia].cname
			glist[gno]={num=0,sst=sst,from=ia,to=ia, cname=cname }
		end
		atomxyz[ia].gno=gno
		glist[gno].num=glist[gno].num+1
		glist[gno].to=ia
	end
	inchain=""
	for ip=1,table.getn(glist) do
		if(inchain~=glist[gno].cname) then
			glist[gno].term="n"
			inchain=glist[gno].cname
		end
	end

	  --- make groups of residues
	groups={}
	sstrec=nil  --- the current group 
	inchain=nil
	for ip=1,table.getn(glist) do
		ginf=glist[ip]
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

	--- for ig=1,table.getn(groups) do print(ig,groups[ig].term) end

	  --- this parameter control how many bones with sizes of
	  --- residues for making groups

	  --- large group threshold 
MAXLEN=15
	  --- group size limit , small groups are leaved out
MINLEN=6
	 --- remove small groups
	 --- except N,C,terminus groups

	 --- first, dissolve small beta-sheet or helix
	 --- and extend loop regions to include them
	for ig=table.getn(groups),1,-1 do
		span=groups[ig]
		len=span[2]-span[1]-1
		if(len<MINLEN) then
	 		print("SMALL",ig,span.name,len,span.sst,span[2]-span[1])
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
					print("marged to",span.name,span.sst)
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
	for ig=table.getn(groups),1,-1 do
		span=groups[ig]
		len=span[2]-span[1]

		  ---	if(span.sst=="H")or(span.sst=="E") then
		if(span.sst~="H") then
			if(len>=MINLEN) then
				if(len>MAXLEN) then
					print("seperate",span[3] or "-",span[1],span[2],len,mid)
					  --- divide two
					mid=math.floor((span[2]+span[1])/2)
					span2={mid+1,span[2],"+"}
					span2.sst="+"
					span[2]=mid
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

function calcavexyz(xyz,pos1,pos2)
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
	return gx,gy,gz
end

function avgvect(p1,p2,atomxyz,allvec)
		local log
		local gx,gy,gz
		gx,gy,gz=calcavexyz(atomxyz,p1,p2)
		grav={gx,gy,gz}
		  --- for all mode
		log=""
		log=log..string.format("%7.3f %7.3f %7.3f ",gx,gy,gz)
		gvec={}
		for im=1,table.getn(allvec) do
			gx,gy,gz=calcavexyz(allvec[im],p1,p2)
			table.insert(gvec,{gx,gy,gz})
		end
		log=log..string.format("%7.3f %7.3f %7.3f ",gx,gy,gz)
		---gx,gy,gz=calcavexyz(avec,1,table.getn(avec))
		---gvec={gx,gy,gz}
		return grav,gvec,log
end
	  --- calc gravity and vector for groups
	bonelist={}
	for ig=1,table.getn(groups) do
		span=groups[ig]

		span.grav,span.gvec,log
			=avgvect(span[1],span[2],atomxyz,allvec)

		  --- also add 1st and last one

		  --- for helix or sheet, apply long range bone

		  --- the control point mid
		mid=math.floor((span[2]+span[1])/2)

--test
if(span.sst=="H") then
mid=span[1]+5
else
mid=span[1]+4
end
		span.mid,span.midvec,log
			=avgvect(mid,mid,atomxyz,allvec)
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
			isp2=mid-1
			iep1=mid+1
		end

		span.ts,span.tsvec
			=avgvect(isp1,isp2,atomxyz,allvec)
		span.te,span.tevec
			=avgvect(iep1,iep2,atomxyz,allvec)
			
---		print(string.format("%s\t%10.4f %10.4f %10.4f %s"
---			,"1st",span.ts[1] ,span.ts[2] ,span.ts[3],span[3] or ""))
----		print(ig,log,isp1,isp2,iep1,iep2,span.name)
---		print(string.format("%s\t%10.4f %10.4f %10.4f %s"
---			,"end",span.te[1] ,span.te[2] ,span.te[3],""))
----		print(span.name,span[1],span[2],span.sst)
	end

	group2bvh(groups,freqlist,molname)
end
-------------

function readNMfile(flcatm,nmode,count)

	local name=string.gsub(flcatm,"#",nmode)
	 --- print(name)
	local fp=io.open(name,"r")
	if( not fp ) then print("xxx no file",name) return end

	line=fp:read()    --- 1 coment 
	line=fp:read()    --- 2 coment 
	line=fp:read()    --- n atoms
	natom=tonumber(string.sub(line,16,-1)) 
	line=fp:read()    --- 
	line=fp:read()    --- freq (1/cm)
	nmode=tonumber(string.sub(line,2,7)) 
	freq=tonumber(string.sub(line,28,37)) 


	line=fp:read()    --- 
	line=fp:read()    --- Atomic ...
	line=fp:read()    --- x  y 

	vectxyz={}

	if not count then
		count=0
	end
	if(count>0) then
		for line in fp:lines() do
			aname=string.sub(line,9,11) 
			if(aname=="CA ") then
				vx=tonumber(string.sub(line,54,61)) 
				vy=tonumber(string.sub(line,62,69)) 
				vz=tonumber(string.sub(line,70,77)) 
				table.insert(vectxyz,{vx,vy,vz})
			end
		end
		print("mode ",nmode,freq)
	else
		print(" ... reading atom info")
		atomxyz={}
		for line in fp:lines() do

			aname=string.sub(line,9,11) 
			if(aname=="CA ") then
				x=tonumber(string.sub(line,27,34)) 
				y=tonumber(string.sub(line,35,42)) 
				z=tonumber(string.sub(line,43,50)) 
				vx=tonumber(string.sub(line,54,61)) 
				vy=tonumber(string.sub(line,62,69)) 
				vz=tonumber(string.sub(line,70,77)) 
				---print(aname,cname,rno,sst,rname,x,y,z,vx,vy,vz)
				arec={x,y,z}
				arec.aname=aname
				arec.ano=tonumber(string.sub(line,1,8))
				arec.cname=string.sub(line,17,17) 
				arec.rno=tonumber(string.sub(line,18,21)) 
				arec.sst=string.sub(line,24,24)
				arec.rname=string.sub(line,12,15) 

				table.insert(atomxyz,arec)
				table.insert(vectxyz,{vx,vy,vz})

			end
			---if( sst == "H") then break end
		end
		print(natom,"atoms",freq,nmode)
	end
	vectxyz.freq=freq
	vectxyz.nmode=nmode

	return vectxyz,atomxyz
end
-----------------
comment=[[
 3162   C   LYS A 402  E    22.714  36.844  13.130     -0.182  -0.051  -0.033   0.192
1---+----|----+----|----+----|----+----|----+----|----+----|----+----|--xx+xxxxx

]]
-----------------


--- for blender, the size of the parent bone
--- the positoin is automatically calculated by blender
--- as the center of all bones

BONELOCx=2
BONELOCy=0
BONELOCz=0

function BoneLoc(x,y,z)
  x=x+BONELOCx
  y=y+BONELOCy
  z=z+BONELOCz
  return x,y,z
end

function group2bvh(groups,freqlist,molname)

	nmode=table.getn(freqlist)
	bonedata={
	nmodefreq=freqlist
	}

	------------------------------

	imode=1

	for ig=1,table.getn(groups) do
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
		mx= (lastlist.xyz[2]+firstlist.xyz[2])/2
		mx= (lastlist.xyz[3]+firstlist.xyz[3])/2
		
			--- for debug purpose
			--- print("F",lastlist[4],firstlist[4],x,y,z,mx,my,mz)

		dx= lastlist.xyz[1]-firstlist.xyz[1]
		dy= lastlist.xyz[2]-firstlist.xyz[2]
		dz= lastlist.xyz[3]-firstlist.xyz[3]
		--	boneone={pos={x=tab[1],y=tab[2],z=tab[3]} }
		boneone={
				axis={x=dx,y=dy,z=dz},
				ctr={x=mx,y=my,z=mz},
		}

		--- the starting pos
		firstp={x=firstlist.xyz[1],y=firstlist.xyz[2],z= firstlist.xyz[3]}
		---midp	={x=mx,	y=	my,z=	 mz}
		lastp={x=lastlist.xyz[1],y=lastlist.xyz[2],z= lastlist.xyz[3]}
		
		boneone.span=span --- from/to data
		boneone.pt1=firstp
		boneone.pt2=lastp
		boneone.pt3={x=x,y=y,z=z} --- midp
		
		tab=midlist
		--- store motion data /translation and rotation
		local vlist={}
		local firstv2,midv2,lastv2
		
		for ig=1,nmode do
			it=ig*3+4	--- the start
			--	vlist[ig]={tab[it],tab[it+1],tab[it+2]}

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
		-------------------
		boneone.name="sm-"..name
		boneone.nmodevect=vlist
		table.insert(bonedata,boneone)
	end
	--------------

		--- calc center of all bone
	nbones=#bonedata
	x=0
	y=0
	z=0
	for ib=1,nbones do
		x=x+bonedata[ib].pt3.x --- pos.x	--TEST
		y=y+bonedata[ib].pt3.y --- pos.y
		z=z+bonedata[ib].pt3.z --- pos.z
	end

	BoneCx=x/nbones
	BoneCy=y/nbones
	BoneCz=z/nbones

	BoneCx,BoneCy,BoneCz=BoneLoc(BoneCx,BoneCy,BoneCz)

	bonedata.pos={x=BoneCx,y=BoneCy,z=BoneCz}

	print("center pos",BoneCx,BoneCy,BoneCz)

		--- mode mixing
	phase={}
	amplitude={}

	for im=1,nmode do
		phase[im]=0+im*0.0001
		bonedata.nmodefreq[im]=freqlist[im]
		---print("f",bonedata.nmodefreq[im],im)
	end

	ntraj=100
	
	amplitude=5	--- 2
	deltat=0.02

	for ib=1,nbones do
		bonedata[ib].motions={}
	end
	
	------ XYZ data util
function XYZadd(a,b,factor)
	if(factor) then
		return {x=a.x+b.x*factor, y=a.y+b.y*factor,	z=a.z+b.z*factor}
	end
	return {x=a.x+b.x, y=a.y+b.y,	z=a.z+b.z}
end

	local sx,sy,sz,mx,my,mz,ex,ey,ez
	local am,mb,sb,mb,eb
	
	for it=1,ntraj do
		t=it*deltat
		for ib=1,nbones do

			axis	=bonedata[ib].axis
			midp	=bonedata[ib].pt3
			firstp=bonedata[ib].pt1
			lastp =bonedata[ib].pt2
			
			mxadd=YotM3RotXyz(0,0,0)
			
			 --- calc 3 division displacement vectors
			 --- make composit vectors for the bone
			mx=0
			my=0
			mz=0
			sx=0
			sy=0
			sz=0
			ex=0
			ey=0
			ez=0
			vects=bonedata[ib].nmodevect
			for im=1,nmode do
				am=amplitude*math.sin(bonedata.nmodefreq[im]*t+phase[im])

				mb=vects[im].vt3
				mx=mx+mb.x*am
				my=my+mb.y*am
				mz=mz+mb.z*am
				fb=vects[im].vt1
				sx=sx+fb.x*am
				sy=sy+fb.y*am
				sz=sz+fb.z*am
				eb=vects[im].vt2
				ex=ex+eb.x*am
				ey=ey+eb.y*am
				ez=ez+eb.z*am
			end

				--- calc the moved location of the 3 divisions 
			firstp2={x=firstp.x+sx,y=firstp.y+sy,z=firstp.z+sz}
			midp2	={x=midp.x+mx,	y=midp.y+my,	z=midp.z+mz}
			lastp2 ={x=lastp.x+ex, y=lastp.y+ey, z=lastp.z+ez}

--- list bones
	if(it==1) then
		len=bonedata[ib].span[2]-bonedata[ib].span[1]+1
		print(ib, len,bonedata[ib].name,bonedata[ib].span[1])
	end
				--- matrix to fit 3 points
			mop,elog=fit3points(firstp,lastp,midp,firstp2,lastp2,midp2)	 
 
			if(elog) then
		print(bonedata[ib].name,it,ib,bonedata[ib].span[1],bonedata[ib].span[2],elog)
		print(string.format("%8.3f %8.3f %8.3f,%8.3f %8.3f %8.3f",sx,sy,sz,mx,my,mz))
		print(string.format("%8.3f %8.3f %8.3f ",ex,ey,ez))
			end 
				if(it==1) and (ib==1) then 
				  -- print(vlist[im][1],im,bonedata.nmodefreq[im])
				end
		  	 --- take the location relative to the center of this object
			tx=mop.center[1]
			ty=mop.center[2]
			tz=mop.center[3]
			
			tx=0
			ty=0
			tz=0

			tx=tx+mop.trn[1]
			ty=ty+mop.trn[2]
			tz=tz+mop.trn[3]
			
			   --- note : BVH format takes ZYX rotation
			rx,ry,rz=YotM3GetZyx(mop.mtx)
		 
	if(elog) then
		print(math.deg(ex),math.deg(ey),math.deg(ez))
	end
		
			  --- use local coordinates of the begining bone
			tx=tx+bonedata[ib].pt1.x	-- -BoneCx
			ty=ty+bonedata[ib].pt1.y	-- -BoneCy
			tz=tz+bonedata[ib].pt1.z	-- -BoneCz

			traj={x=tx,y=ty,z=tz,xrot=rx,yrot=ry,zrot=rz}
			bonedata[ib].motions[it]=traj
		end
	end

	noaxis={x=5,y=0,z=0}

		-----------
		
	nframes=ntraj
	tframe=0.014
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
	for im=1,nframes do
		for ib=1,nbones do
			locrot=bonedata[ib].motions[im]

			io.write(string.format("%.5f %.5f %.5f",locrot.x-BoneCx,locrot.z-BoneCz,-(locrot.y-BoneCy)))
																						--- SWAP ZY
			xrot=locrot.xrot or 0
			yrot=locrot.yrot or 0
			zrot=locrot.zrot or 0
			io.write(string.format(" %.3f %.3f %.3f ",xrot,zrot,-yrot))
		end
		io.write("\n") --- end of skelton
	end
	io.write("\n") --- end of skelton

end

--------------------
main()
