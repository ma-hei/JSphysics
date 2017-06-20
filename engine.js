
var RigidBody = function(){
	
	this.mass = 10;
	this.I_body = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.I_body_inv = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);	
	this.x = [0,0,1];
	this.q = [0,0,0,0];
	this.P = [0,0,0];
	this.L = [0,0,0];
	this.I_inv = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.R = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.v = [0,0,0];
	this.omega = [0,0,0];
	this.force = [0,0,0];
	this.torque = [0,0,0];
}


var Simulation = function(){

	this.n_bodies = 1;
	this.state_size = 18;
	this.rigid_bodies = new Array(this.n_bodies);
	
}

Simulation.prototype.quaterion_to_matrix = function(q){
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

Simulation.prototype.init_states = function(){

	for (var i=0;i<this.n_bodies;i++){
		this.rigid_bodies[i] = new RigidBody();
	}
	
}

Simulation.prototype.state_to_array = function(rb, y, idx){

	y[idx+0] = rb.x[0];
	y[idx+1] = rb.x[1];
	y[idx+2] = rb.x[2];

	y[idx+3] = rb.q.r;
	y[idx+4] = rb.q.i;
	y[idx+5] = rb.q.j;
	y[idx+6] = rb.q.k;

	y[idx+7] = rb.P[0];
	y[idx+8] = rb.P[1];
	y[idx+9] = rb.P[2];

	y[idx+10] = rb.L[0];
	y[idx+11] = rb.L[1];
	y[idx+12] = rb.L[2];
}

Simulation.prototype.array_to_state = function(rb, y, idx){

	rb.x[0] = y[idx+0];
	rb.x[1] = y[idx+1];
	rb.x[2] = y[idx+2];

	rb.q.r = y[idx+3];
	rb.q.i = y[idx+4];
	rb.q.j = y[idx+5];
	rb.q.k = y[idx+6];

	rb.P[0] = y[idx+7];
	rb.P[1] = y[idx+8];
	rb.P[2] = y[idx+9];
	
	rb.L[0] = y[idx+10];
	rb.L[1] = y[idx+11];
	rb.L[2] = y[idx+12];

	rb.v[0] = rb.P[0]/rb.mass;
	rb.v[1] = rb.P[1]/rb.mass;
	rb.v[2] = rb.P[2]/rb.mass;

	rb.R = this.quaterion_to_matrix(rb.q);
	
	rb.I_inv = math.multiply(math.multiply(rb.R,rb.I_body_inv), math.transpose(rb.R));

	rb.omega = math.multiply(rb.I_inv,rb.L);

}

Simulation.prototype.bodies_to_array = function(x){

	for (var i=0; i<this.n_bodies; i++){
		this.state_to_array(this.rigid_bodies[i], x, i*this.state_size);
	}
}

Simulation.prototype.array_to_bodies = function(x){
	for (var i=0;i<this.n_bodies; i++){
		this.array_to_state(this.rigid_bodies[i], x, i*this.state_size);
	}
}

Simulation.prototype.quaternion_product(a,b){
	//[s1,v1] * [s2,v2] = s1s2 - v1.v2, s1v2+s2v1 + v1xv2]
	
}

Simulation.prototype.DdtStateToArray(rb, xdot, idx){
	xdot[idx+0] = rb.v[0];
	xdot[idx+1] = rb.v[1];
	xdot[idx+2] = rb.v[2];

	
}

Simulation.prototype.Dxdt(t, x, xdot){

	this.array_to_bodies(x);
	for (var i=0; i<this.n_bodies; i++){
		this.compute_force_and_torque(t, this.rigid_bodies[i]);
	}

}

Simulation.prototype.ode(x0, xFinal){
}

Simulation.prototype.run_simulation = function(){

	var x0 = new Array(this.n_bodies * this.state_size);
	var xFinal = new Array(this.n_bodies * this.state_size);

	this.init_states();
	this.bodies_to_array(xFinal);
	for (var t=0; t<10.0; t+=1/24){
		for (var i=0;i<this.state_size*this.n_bodies;i++){
			x0[i] = xFinal[i];
		}
		//this.ode(x0, xFinal, this.state_size*this.n_bodies, t, t+1/
		this.ode(x0, xFinal);
		this.array_to_bodies(xFinal);
	}
}

