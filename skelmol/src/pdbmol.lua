
  --- lua5.0 5.3 compatible
  --- pdb file read 
  --- ver 1.01 multiple model support 2017 April

if(not msCreateMol) then
	msCreateMol=function()
		pdb={
			sequence={}
		}
		return pdb
	end
	msCreateAtom=function(mol,atom)
		local seq
		table.insert(mol,atom)
		seq=mol.sequence[atom.cname]
		if(not seq) then
			seq={}
			mol.sequence[atom.cname]=seq
		end
		seq[atom.resno]=atom.resname
		if(not atom.bond) then
			atom.bond={}
		end
		return atom
	end

end

  --- test if two atoms are bonded
msIsBond=function(atom1,atom2)
	local bond=atom1.bond
	if(not bond) then
		return
	end
	for idx,cnt in pairs(bond) do
		if(cnt==atom2) then
			return atom2
		end
	end
end

msSetBondTable=function(atom1,atom2)
	local bond=atom1.bond
	if(not bond) then
		bond={}
		atom1.bond=bond
	end
	local bond2=atom2.bond
	if(not bond2) then
		bond2={}
		atom2.bond=bond2
	end
	if(msIsBond(atom1,atom2)) then
		return
	end
	table.insert(bond,atom2)
	table.insert(bond2,atom1)
	return #bond
end

  --- make a clone copy, sharing contents
msCloneMol=function(mol)
	local one
	local newmol={}

	for lp=1,getn(mol) do
		one={}
		for idx,cnt in pairs(mol[lp]) do
			one[idx]=cnt
		end
		newmol[lp]=one
	end

	for idx,cnt in pairs(mol) do
		if(type(idx)=="string") then
			newmol[idx]=cnt
		end
	end
	return newmol
end


function msReadPdb(file)
	local pdb,aname
	local resno,rname
	local cname

	if(file) then
		if(not io.input(file)) then
			error("pdb file "..file)
		end
	end

	pdb=msCreateMol()
	line=""
	while(line) do
		line=io.read()
		if(not line) then break end
		if( string.find(line,"^ATOM")
		 or string.find(line,"^HETATM")) then
			cname=string.sub(line,22,22)	
			rname=string.sub(line,18,20)	

			atom={}
			atom.pdb=tonumber(string.sub(line, 7,12))	
			atom.aname=string.sub(line,13,16)	   ---x
			atom.resno=tonumber(string.sub(line,23,26))	
			if(not atom.resno) then
				print(" XXX error in residue number ")
				print(string.sub(line,23,26))	
			end
			atom.resname=rname
			atom.cname=cname
			atom.x=tonumber(string.sub(line,31,38))	
			atom.y=tonumber(string.sub(line,39,46))	
			atom.z=tonumber(string.sub(line,47,54))	
			atom.occu=tonumber(string.sub(line,55,60))	
			atom.temp=tonumber(string.sub(line,61,65))	
			msCreateAtom(pdb,atom)
			  --- support multi model file
		elseif( string.find(line,"^ENDMDL")) then
			break
		elseif( string.find(line,"^MODEL")) then
			pdb.modelno=tonumber(string.sub(line,9,15))
		end
	end
	return pdb
end

comment=[[#ifdef  COMMENT
ATOM      1  CA  GLY A   1  0   14.021  46.337  31.119  1.00 34.55   1
%%%%%%      %%%% %%%  ####    +###.###        +###.###      +##.##          %%
      #####     %    %    %           +###.###        +##.##       ###  %%%%  %%
1---+----|----+----|----+----|----+----|----+----|----+----|----+----|--xx+xxxxx

         1..(6) "ATOM"
         7..(5) #number                 31..(8) #x
        l3..(4) name                    39..(8) #y
        l7      alternate location      47..(8) #z
        18..(3) residue name            55..(6) #occupancy
        22      chain name              61..(6) #temperature factor
        23..(4) #residue number         68..(3) #footnote number
        27      res. insertion code
       x73..(4) segID   (new)
       x77..(2) element (new)
       x79..(2) charged (new)
#endif
]]


------------------------------------

function msWritePdb(mol)
	local seq
	for ia=1,table.getn(mol) do
		atom=mol[ia]
		name=atom.aname or "XX__"
		if(string.sub(name,1,1)=="_") then
			name=" "..string.sub(name,2)
		end

		seq=mol.sequence[atom.cname]
		if(seq) then
			res=seq[atom.resno] or "---"
		else
			res="---"
		end
		chain=atom.cname or ""
		resno=atom.resno or 1
		io.write(string.format(
     "ATOM %6d %4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f%6.2f",
			atom.pdb or ia,
			name,res,chain,resno,
			atom.x,atom.y,atom.z,atom.occu or 1,atom.temp or 0 ))
		io.write("        \n")
	end
end

-----------------------

function msPdbBond(mol)
	local ires,bmax2,bmax
	local dx,dy,dz

	bmax=1.9
	bmax2=bmax*bmax

	nmax=table.maxn(mol)
	for idx=1,nmax do
		one=mol[idx]

		--print("dd",one.pdb,one.res)
		ires=one.resno
		for idx2=idx+1,nmax do
			two=mol[idx2]
			if ires and (two.resno>ires+1) then
				  --- print("NN",two.res,ires)
				break
			end
				  --- print("cc",two.pdb,two.res,ires)
			dx=one.x-two.x
			dy=one.y-two.y
			dz=one.z-two.z
			dx=dx*dx+dy*dy+dz*dz
			if(dx<bmax2) then
				msSetBondTable(one,two)
			end
		end
	end
end

------------------------------------
