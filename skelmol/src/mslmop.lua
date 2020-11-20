
---
--- Mop  --- twist bug fixed 2013 Nov

print("--- in Mop")

function compat()
	print("XXX pease avoid old func sin/cos")
	sin=function(xx) return math.sin(math.rad(xx)) end
	cos=function(xx) return math.cos(math.rad(xx)) end
	acos=function(xx) return math.deg(math.acos(xx)) end
	atan2=function(xx,yy) return math.deg(math.atan2(xx,yy)) end
end
function oldsin(x) compat() return sin(x) end
function oldcos(x) compat() return cos(x) end
sin=oldsin
cos=oldcos

Mop={}

function Mop.new()
	local mop={}

	mop.trn={0,0,0}
	mop.mtx={1,0,0, 0,1,0,  0,0,1}
	  --- no matrix by default
	Mop.proc.setup(mop)
	return mop
end

Mop.proc={}

Mop.proc.setup=function(mop)
	for idx,cnt in pairs(Mop.proc) do
		mop[idx]=cnt
	end
end

	   --- rotation -x-z-x
Mop.proc.euler=function(mop,a,b,c)
	local cr,ct,cx,cy,cz
	local s1,s2,s3,c1,c2,c3
	local mtx11,mtx12,mtx13,mtx21,mtx22,mtx23,mtx31,mtx32,mtx33
	local tx,ty,tz

	s1=sin(a)
	s2=sin(b)
	s3=sin(c)
	c1=cos(a)
	c2=cos(b)
	c3=cos(c)

	mtx11=-s1*c2*s3+c1*c3
	mtx12=c1*c2*s3+s1*c3
	mtx13=s2*s3
	mtx21=-s1*c2*c3-c1*s3
	mtx22=c1*c2*c3-s1*s3
	mtx23=s2*c3
	mtx31=s1*s2
	mtx32=-c1*s2
	mtx33=c2

	mop.mtx={mtx11,mtx12,mtx13,mtx21,mtx22,mtx23,mtx31,mtx32,mtx33}
end
Mop.proc.twist=function(mop,vx,vy,vz,t)
	local dn
	local c3,s3,c3c

	dn=math.sqrt(vx*vx+vy*vy+vz*vz)
	if (dn<1e-7) then
		print("too little vector",vx,vy,vz,t)
		mop:twist(vx+1e-7,vy,vz,t)
		print("mop",mop.mtx[1])
		return mop
	end
	vx=vx/dn
	vy=vy/dn
	vz=vz/dn
	c3=cos(t)
	s3=sin(t)
	c3c=1.0-c3
	mtx11 = vx*vx*c3c  + c3 
	mtx12 = vx*vy*c3c  + vz*s3
	mtx13 = vx*vz*c3c  - vy*s3 
	mtx21 = vy*vx*c3c  - vz*s3 
	mtx22 = vy*vy*c3c  + c3   
	mtx23 = vy*vz*c3c  + vx*s3 
	mtx31 = vz*vx*c3c  + vy*s3
	mtx32 = vz*vy*c3c  - vx*s3
	mtx33 = vz*vz*c3c  + c3  

	mop.mtx={mtx11,mtx12,mtx13,mtx21,mtx22,mtx23,mtx31,mtx32,mtx33}
end
Mop.proc.calc=function(mop,x,y,z) 
	local cx,cy,cz,tx,ty,tz
	local mtx

	if(mop.center) then
		cx=mop.center[1]
		cy=mop.center[2]
		cz=mop.center[3]
	else
		cx=0
		cy=0
		cz=0
	end

	mtx=mop.mtx
	if(mtx) then
		tx=x-cx
		ty=y-cy
		tz=z-cz
		x=cx+mtx[1]*tx+mtx[2]*ty+mtx[3]*tz
		y=cy+mtx[4]*tx+mtx[5]*ty+mtx[6]*tz
		z=cz+mtx[7]*tx+mtx[8]*ty+mtx[9]*tz
	end
	if(mop.trn) then
		x=x+mop.trn[1]
		y=y+mop.trn[2]
		z=z+mop.trn[3]
	end
	return x,y,z
end

Mop.proc.setcam=function(mop) 
	PufCamMtx_pier(mop.mtx[1],mop.mtx[2],mop.mtx[3],
			mop.mtx[4],mop.mtx[5],mop.mtx[6],
			mop.mtx[7],mop.mtx[8],mop.mtx[9],
			mop.trn[1],mop.trn[2],mop.trn[3])
end
Mop.proc.unsetcam=function()
	PufCamMtx_pier()
end

Mop.proc.calcxyz=function(mop,t) 
	t.x,t.y,t.z=mop:calc(t.x,t.y,t.z)
end

	--------------------

  --- from cdyel.lua
  --- euler angle

function deg360(deg) 
	if(deg>=0) then
		deg=mod(deg,360)
	else
		deg=360+mod(deg,360)
	end
	return deg
end


----- ** bug fixed 1999 Aug 10 y.ueno **/
----- ** z-x-z type euler angle is used here	**/


function ComlinRange(a1,a2)
	a1=deg360(a1)
	a2=deg360(a2)
	if(a1>180) then
		a1=a1-180
		a2=a2-180
		if(a2<0) then
			a2=a2+360
		end
	end
	return a1,a2
end

  --- ** note there are 2 possible solution for euler angle **/
  --- ** due to handedness **/
  --- ** it retuns 0 if two matix is identical or equal to unity **/
  --- ** ... not yet **/

function Colin2Euler(cm6arg)
	local	angA,angB,angC,tilta,tiltb;
	local	cosA,cosB,cosC,sinA,sinB,sinC;
	local	ang12,ang21,ang13,ang31,ang23,ang32;
	local	cm6test,cm6;
	local idx,idis;
	local	testfloat;
	local ea2,ea3
	  --- ** very important to use float to reduce precisions **/

	  --- ** first clear matrices **/
	ea2={}
	ea3={}

	cm6={}
	cm6[1],cm6[2]=ComlinRange(cm6arg[1],cm6arg[2])
	cm6[3],cm6[4]=ComlinRange(cm6arg[3],cm6arg[4])
	cm6[5],cm6[6]=ComlinRange(cm6arg[5],cm6arg[6])
	

	ang12=cm6[1];
	ang21=cm6[2];

	ang13=cm6[3];
	ang31=cm6[4];

	ang23=cm6[5];
	ang32=cm6[6];

	angA=ang13-ang12;
	angB=ang23-ang21;
	angC=ang32-ang31;


	  --- ** we can try to avoid case with angA==0 or angB==0 **/
	if(angA==0) then
		 --- ** print(" zero sin A \n"); /**/
		ang13=ang13+0.1;
		angA=0.1;
	end
	if(angB==0) then
		 --- ** fprintf(stderr," zero sin B \n"); /**/
		ang23=ang23+0.1;
		angA=0.1;
	end

	cosA=cos(angA)
	cosB=cos(angB)
	cosC=cos(angC)
	sinA=sin(angA)
	sinB=sin(angB)
	sinC=sin(angC)

	testfloat=sinA*sinB
	if(testfloat==0) then
		print(format(" zero sin %f %f \n",angA,angB))
		return
	end

	testfloat=((cosC-cosA*cosB)/(sinA*sinB));
	if((abs(testfloat)>1.) ) then
		print(format(" funny anges for acos() %f \n",testfloat))
		return
	end

	tilta=acos(testfloat);

	  --- ** fprintf(stderr," angles %6.2f %6.2f %6.2f", ang12, ang21,ang31); **/
	  --- ** fprintf(stderr," a=%7.2f b=%7.2f c=%7.2f\n",angA,angB,angC); **/


	testfloat= sinB/sinC*sin(tilta);
	if((abs(testfloat)>1.) ) then
		print(format(" funny anges for asin() %f \n",testfloat))
		return
	end
	tiltb=asin(testfloat);

	  --- ** fprintf(stderr," calculated plane tilt anges %6.2f %6.2f\n"
	  ---	,tilta,tiltb);
	  --- ** tiltb angle may be calculated by another	**/
	  --- ** formula but has no sign information		**/
	  --- **	test=((cosC*cosA-cosB)/(sinA*sinC));
	  --- **	test=RAD2DEG*acos(test);
	  --- **	fprintf(stderr," another value  %6.2f %6.2f\n",tiltb,test);
	  --- ** **/

	  --- ** we prefer beta angle in 0..180 **/
	  --- ** so flip alpha and gamma	**/
	if(tiltb<0) then
		tiltb= -tiltb;
		ang13=ang13+180;
		ang31=ang31+180;
	end

	ea2[1]= deg360( ang12);
	ea2[2]= deg360(tilta);
	ea2[3]= deg360(-ang21);
	ea3[1]= deg360( ang13);

	ea3[2]= deg360(tiltb);
	ea3[3]= deg360(-ang31);

	  --- ** then test reversal to common line **/
	  --- ** to find impossible angles 	**/

DEGDELTA = 20

	cm6test=Euler2Colin(ea2,ea3)
	if(not cm6test) then
 		print(format(" not a unique euler"))
	else
		---- for idx=1,6 do
		idx=1 while(idx<=6) do
			idis=abs(cm6test[idx]-cm6[idx]);
			if(idis>DEGDELTA) then
			print(format(" not a uniqe cm %f %f \n",cm6test[idx],cm6[idx]))
			end
		idx=idx+1 end
	end

	  --- ** check if two beta angles are same or 0 **/
	idis=abs(ea2[2])
	if(idis<1) then
 		printf(format(" zero tilt2 %f \n",ea2[1]))
	end
	idis=abs(ea3[2])
	if(idis<1) then
 		printf(format(" zero tilt3 %f \n",ea3[1]))
	end

	return ea2,ea3
end

function Euler2Colin(ea2,ea3)
	local aone,atwo,stat
	local mea1,mea2,mea3,cm6
	mea1=EulerMtx(0,0,0)
	mea2=EulerMtx(ea2[1],ea2[2],ea2[3])
	mea3=EulerMtx(ea3[1],ea3[2],ea3[3])

	stat=1
	aone,atwo=MtxColin(mea1,mea2)
	if(not aone) then
		stat=nil
	end

	cm6={}
	cm6[1]=aone
	cm6[2]=atwo

	aone,atwo=MtxColin(mea1,mea3)
	if(not aone) then
		stat=nil
	end
	cm6[3]=aone
	cm6[4]=atwo

	aone,atwo=MtxColin(mea2,mea3)
	if(not aone) then
		stat=nil
	end
	cm6[5]=aone
	cm6[6]=atwo

	return  cm6

end



---
--- Mop 

		--- Mop.proc.euler=function(mop,a,b,c)
	EulerMtx=function(a,b,c)
		local cr,ct,cx,cy,cz
		local s1,s2,s3,c1,c2,c3
		local mtx11,mtx12,mtx13,mtx21,mtx22,mtx23,mtx31,mtx32,mtx33
		local tx,ty,tz

		s1=sin(a)
		s2=sin(b)
		s3=sin(c)
		c1=cos(a)
		c2=cos(b)
		c3=cos(c)

		mtx11=-s1*c2*s3+c1*c3
		mtx12=c1*c2*s3+s1*c3
		mtx13=s2*s3
		mtx21=-s1*c2*c3-c1*s3
		mtx22=c1*c2*c3-s1*s3
		mtx23=s2*c3
		mtx31=s1*s2
		mtx32=-c1*s2
		mtx33=c2

		return {mtx11,mtx12,mtx13,mtx21,mtx22,mtx23,mtx31,mtx32,mtx33}
	end
	TwistMtx_OLD=function(dx,dy,dz,t)
		local vx,vy,vz,dn
		local c3,s3,c3c

		dn=sqrt(vx*vx+vy*vy+vz*vz)
		if (dn<1e-7) then
			print("too little vector")
			return
		end
		vx=vx/dn
		vy=vy/dn
		vz=vz/dn
		c3=cos(t3)
		s3=sin(t3)
		c3c=1.0-C3
		mtx11 = vx*vx*c3c  + c3 
		mtx12 = vx*vy*c3c  + vz*s3
		mtx13 = vx*vz*c3c  - vy*s3 
		mtx21 = vy*vx*c3c  - vz*s3 
		mtx22 = vy*vy*c3c  + c3   
		mtx23 = vy*vz*c3c  + vx*s3 
		mtx31 = vz*vx*c3c  + vy*s3
		mtx32 = vz*vy*c3c  - vx*s3
		mtx33 = vz*vz*c3c  + c3  

		return {mtx11,mtx12,mtx13,mtx21,mtx22,mtx23,mtx31,mtx32,mtx33}
	end


  --- ** acos() without error, always return valid value **/
function acosWOE(val)
	if(val<= -1) then
		return 180
	elseif(val>=1.) then
		return 0
	end
	return acos(val)
end


  --- ** get an euler angle from a matrix **/
function MtxEuler(dest)
	local	a,b,c,s2,c2;

					-- m[9]=cosB
	b=acosWOE(dest[9])
	s2=sin(b)	--- ** but not correct b **/

	  --- ** if beta is vey small, it probably cause error **/
	if(b<1e-5) then
		s2=0;
	end

	if(s2==0) then
		  --- ** this case, only (a+c) value is found **/
		  --- ** so let c=0 **/
		b=0
		c=0
					---M11 m[1]= -sinA*cosB*0+cosA*1
		a=acosWOE(dest[1])
					---M12 m[2]= cosA*cosB*0+sinA*1
		if((dest[2])<0) then	--- * asin() unchanges +-sign **/
			a=360-a;
		end
	else 
					---M32 m[8]= -cosA*sinB
		a=acosWOE(-dest[8]/s2)
				--- ** value is between 0..PI **/
					---M31 m[7]=sinA*sinB
		if((dest[7]/s2)<0) then --- * asin() unchanges +-sign **/
			a=360-a;
		end
					---M23 m[6]= sinB*cosC
		c=acosWOE(dest[6]/s2)
					---M13 m[3]= sinB*sinC
		if((dest[3]/s2)<0) then
			c=360-c
		end
		if(s2<0) then
			b=360-b
		end
	end
	a=deg360(a)
	b=deg360(b)
	c=deg360(c)

	if(a<0) then
		print("NEG",a)
	end
	return a,b,c
end

function YotM3Trans(m)
	return {m[1],m[4],m[7],
		m[2],m[5],m[8],
		m[3],m[6],m[9]}
end

function YotM3xV3(mtx,v)
	local x,y,z
	x=mtx[1]*v[1]+mtx[2]*v[2]+mtx[3]*v[3]
	y=mtx[4]*v[1]+mtx[5]*v[2]+mtx[6]*v[3]
	z=mtx[7]*v[1]+mtx[8]*v[2]+mtx[9]*v[3]
	return {x,y,z}
end

function YotM3xM3(ma,mb)
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


function MtxColin(mtx1,mtx2)
	local n1,n2,v1,v2,cm12,zaxis
	local a1,a2,dlen
	local invMtx
	local stat

	stat=TRUE;
	  --- ** first calculate 2 plane vector **/
	  --- ** get the world coordinates from a local coordinate of the plane **/
	zaxis={0,0,1}
	invMtx=YotM3Trans(mtx1)
	n1=YotM3xV3(invMtx,zaxis)
	invMtx=YotM3Trans(mtx2)
	n2=YotM3xV3(invMtx,zaxis)

	cm12={}
	cm12[1]= n1[2]*n2[3] - n1[3]*n2[2]
	cm12[2]= n1[3]*n2[1] - n1[1]*n2[3]
	cm12[3]= n1[1]*n2[2] - n1[2]*n2[1]

	dlen=cm12[1]*cm12[1]+cm12[2]*cm12[2]+cm12[3]*cm12[3]
	dlen=sqrt(dlen)
	if( dlen==0 ) then
		stat=FALSE;
		  --- ** it can be any vector except z axis **/
		cm12[1]=n1[2]
		cm12[2]=n1[3]
		cm12[3]=n1[1]
	else
		cm12[1]=cm12[1]/dlen
		cm12[2]=cm12[2]/dlen
		cm12[3]=cm12[3]/dlen
	end

	v1=YotM3xV3(mtx1,cm12);
	v2=YotM3xV3(mtx2,cm12);

	a1= atan2(v1[2],v1[1])
	if(a1<0) then
		a1=a1+360
	end
	a2= atan2(v2[2],v2[1])
	if(a2<0) then
		a2=a2+360
	end

		--- ** range of the angle one is [0..PI] **/
		--- ** so, take reverse one (conjugate) **/
	if(a1>=180) then
		a1=a1-180
		a2=a2-180
		if(a2<0) then
			a2=a2+360
		end
	end

	  --- ** a1= acos(v1[1]); **/
	  --- ** if(asin(v1[2])<0) a1=PI*2-a1; **/

	return a1,a2
end

