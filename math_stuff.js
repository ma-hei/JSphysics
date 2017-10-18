function quaterion_to_matrix(q){
 	//q = s, vx, vy, vz
	var s = q[0];
	var vx = q[1];
	var vy = q[2];
	var vz = q[3];

	var m11 = 1-2*vy*vy-2*vz*vz;
	var m12 = 2*vx*vy-2*s*vz;
	var m13 = 2*vx*vz+2*s*vy;

	var m21 = 2*vx*vy+2*s*vz;
	var m22 = 1-2*vx*vx-2*vz*vz;
	var m23 = 2*vy*vz-2*s*vx;

	var m31 = 2*vx*vz-2*s*vy;
	var m32 = 2*vy*vz+2*s*vx;
	var m33 = 1-2*vx*vx-2*vy*vy;
	var m = math.matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]);

	return m;
 
}
