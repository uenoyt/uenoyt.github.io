
--------------------------------------------------------------
--- fiting 2 vectors
----

dofile("mslmop.lua")

function YotM3RotXyz(xrot,yrot,zrot)
  local cosx,sinx,cosy,siny,cosz,sinz
  local mxX,mxY,mxZ,myX,myY,myZ,mzX,mzY,mzZ

  xrot=math.rad(xrot)
  yrot=math.rad(yrot)
  zrot=math.rad(zrot)
	cosx=math.cos(xrot)    sinx=math.sin(xrot)
	cosy=math.cos(yrot)    siny=math.sin(yrot)
	cosz=math.cos(zrot)    sinz=math.sin(zrot)
		--/* make x,y,z rotation matrix */
	mxX= cosy*cosz;
	mxY= sinx*siny*cosz-cosx*sinz;
	mxZ= cosx*siny*cosz+sinx*sinz;
	myX= cosy*sinz;
	myY= sinx*siny*sinz+cosx*cosz;
	myZ= cosx*siny*sinz-sinx*cosz;
	mzX= -siny;
	mzY= cosy*sinx;
	mzZ= cosx*cosy;
  return {mxX,mxY,mxZ,myX,myY,myZ,mzX,mzY,mzZ}
end

---- calc xyz rotation angles from a matrix
function YotM3GetXyz(mtx)
  local 	ry,rx,rz,sinz,cosz,siny,cosy,sinx,cosx
  local   mxX,mxY,mxZ,myX,myY,myZ,mzX,mzY,mzZ
    = mtx[1],mtx[2],mtx[3],mtx[4],mtx[5],mtx[6],mtx[7],mtx[8],mtx[9]

  local PI=math.pi
		--/** rr.zY= cosy*sinx; **/
		--/** rr.zZ= cosx*cosy; **/
	rx=math.atan2(mzY,mzZ)	  --/** x= -90..90 **/
	if(rx>(PI/2)) then
    rx=rx-PI
	elseif(rx<(-PI/2)) then
    rx=rx+PI
  end
				  --/** z= -180..180 **/
	sinx=math.sin(rx)
	cosx=math.cos(rx)

	cosz=cosx*myY-sinx*myZ
	sinz=-cosx*mxY+sinx*mxZ
	rz=math.acos(cosz)
	if(sinz<0) then
    rz=-rz
  end  -- /** y= -180..180 **/
	siny= -mzX
	ry=math.asin(siny)
				--/** cosy= mtx->yX/sinz; **/
				--/** rr.zZ= cosx*cosy; **/
	if(cosx~=0) then
		cosy= mzZ/cosx
	else
		cosy= mzY/sinx
  end
	ry=math.acos(cosy)
	siny= -mzX
	if(siny<0) then
    ry=-ry
  end
  return rx,ry,rz
end

  ---- calc ZYX rotation angles from 3x3 matrix
  
function YotM3GetZyx(m)
  local minv=YotM3Trans(m)
  local irx,iry,irz=YotM3GetXyz(minv)
  irx,iry,irz=math.deg(irx),math.deg(iry),math.deg(irz)
  return -irx,-iry,-irz
end


function fit3points(pa1,pa2,pa3,pb1,pb2,pb3)
  local gx1,gy1,gz1,gx2,gy2,gz2
--- cc	!! use first 3 point
	vx1=pa1.x-pa2.x
	vy1=pa1.y-pa2.y
	vz1=pa1.z-pa2.z
	ux1=pa2.x-pa3.x
	uy1=pa2.y-pa3.y
	uz1=pa2.z-pa3.z

	vx2=pb1.x-pb2.x
	vy2=pb1.y-pb2.y
	vz2=pb1.z-pb2.z
	ux2=pb2.x-pb3.x
	uy2=pb2.y-pb3.y
	uz2=pb2.z-pb3.z

	  --- center of the 3 atoms
	gx1= (pa1.x+pa2.x+pa3.x)/3
	gy1= (pa1.y+pa2.y+pa3.y)/3
	gz1= (pa1.z+pa2.z+pa3.z)/3

	gx2= (pb1.x+pb2.x+pb3.x)/3
	gy2= (pb1.y+pb2.y+pb3.y)/3
	gz2= (pb1.z+pb2.z+pb3.z)/3

	print(string.format(" grav1 %9g %9g %9g ",gx1,gy1,gz1))
	print(string.format(" grav2 %9g %9g %9g ",gx2,gy2,gz2))

	  --!! use first 3 point

	print(string.format("fit vector 1 %9g %9g %9g  %9g %9g %9g "
		,vx1,vy1,vz1,ux1,uy1,uz1))
	print(string.format("fit vector 2 %9g %9g %9g  %9g %9g %9g "
		,vx2,vy2,vz2,ux2,uy2,uz2))

	mop=fit2vec(vx1,vy1,vz1,ux1,uy1,uz1
     			      ,vx2,vy2,vz2,ux2,uy2,uz2)

  if(nil) then  --- if rotation center is origin
    cx,cy,cz=mop:calc(gx1,gy1,gz1)  
    mop.trn[1]=gx2-cx
    mop.trn[2]=gy2-cy
    mop.trn[3]=gz2-cz
  else
    mop.center={gx1,gy1,gz1}
    mop.trn[1]=gx2-gx1
    mop.trn[2]=gy2-gy1
    mop.trn[3]=gz2-gz1
  end
return mop
end


function v3mul(fx,fy,fz,gx,gy,gz)
	local ax,ay,ay
        ax=fy*gz-fz*gy
        ay=fz*gx-fx*gz
        az=fx*gy-fy*gx
	return ax,ay,az
end

function v3ang(fx,fy,fz,gx,gy,gz)
        local    fxr,fyr,fzr, gxr,gyr,gzr
        local    fr2,fr,frr,  gr2,gr,grr, cst,ang

        gr2=gx*gx+gy*gy+gz*gz
        gr=math.sqrt(gr2)
        if(gr < 1e-16) then
          print(" XXX error v3ang 1st vector 0 ")
          gr=1.
        end
        gxr=gx/gr
        gyr=gy/gr
        gzr=gz/gr

        fr2=fx*fx + fy*fy + fz*fz
        fr=math.sqrt(fr2)
        if(fr<1e-16) then
          print(" XXX error v3ang 2nd vector 0 ")
          fr=1.
        end
        fxr=fx/fr
        fyr=fy/fr
        fzr=fz/fr
        cst=fxr*gxr+fyr*gyr+fzr*gzr

        if (math.abs(cst)>0.9999) then
	  cst=math.max(-0.9999,math.min(0.9999,cst))
          print( " XXX warning. angle is almost linear.")
        end

        ang=math.deg(math.acos(cst))
        return ang
end

function m3xm3(ma,mb)
        local mc={
         ma[1]*mb[1]+ma[2]*mb[4]+ma[3]*mb[7],
         ma[1]*mb[2]+ma[2]*mb[5]+ma[3]*mb[8],
         ma[1]*mb[3]+ma[2]*mb[6]+ma[3]*mb[9],

         ma[4]*mb[1]+ma[5]*mb[4]+ma[6]*mb[7],
         ma[4]*mb[2]+ma[5]*mb[5]+ma[6]*mb[8],
         ma[4]*mb[3]+ma[5]*mb[6]+ma[6]*mb[9],

         ma[7]*mb[1]+ma[8]*mb[4]+ma[9]*mb[7],
         ma[7]*mb[2]+ma[8]*mb[5]+ma[9]*mb[8],
         ma[7]*mb[3]+ma[8]*mb[6]+ma[9]*mb[9]
        }
        return mc
end

function fit2vec(vx1,vy1,vz1,ux1,uy1,uz1
     			      ,vx2,vy2,vz2,ux2,uy2,uz2)
	  --!! fit two vector
	  --!!   v1 v2 is a major vector to be fitted
	  --!!   u1 u2 is a sub vector which may not be parpendicular to v1 v2
	  --!! pepair sub vector o1 o2 vector which is normal to v1 v2
	  --!! 
	  --!! fit v1 v2 : first get a v3 vector which is normal to v1 v2
	  --!!           : then rotate around v3  angle v1^v2
	  --!! get new o1  applying the rotation
	  --!! fit o1' o2: roate around v2 angle  v1'^v2
	  --!!	       ( v1 and v2 is normal to o3 )
	local ox1,oy1,oz1, ox2,oy2,oz2,t
	local mtx2,mtx3

	ox1,oy1,oz1= v3mul(vx1,vy1,vz1,ux1,uy1,uz1)
	ox2,oy2,oz2= v3mul(vx2,vy2,vz2,ux2,uy2,uz2)

	  --- print( ' sub vecotr 1 ',ox1,oy1,oz1)
	  --- print( ' sub vecotr 2 ',ox2,oy2,oz2)

	vx3,vy3,vz3= v3mul(vx1,vy1,vz1,vx2,vy2,vz2)
	t=-v3ang(vx1,vy1,vz1,vx2,vy2,vz2)

	  --- print( ' 1st rotation angle ',t,' around vector',vx3,vy3,vz3)
	mtx2=Mop:new()
	mtx2:twist(vx3,vy3,vz3,t)

	ox1,oy1,oz1=mtx2:calc(ox1,oy1,oz1)

	t=-v3ang(ox1,oy1,oz1,ox2,oy2,oz2)

	  --- print(' 2nd rotation angle ',t,' around vector',vx2,vy2,vz2)
	mtx3=Mop:new()
	mtx3:twist(vx2,vy2,vz2,t)

	mtx3.mtx=m3xm3(mtx3.mtx,mtx2.mtx)

	return mtx3
end 

---------------------

