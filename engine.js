
var RigidBody = function(){
	
	this.mass = 10;
	this.I_body = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.I_body_inv = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);	
	this.x = [20,20,0];
	this.q = [0,0,0,1];
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
	this.state_size = 13;
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

	y[idx+3] = rb.q[0];
	y[idx+4] = rb.q[1];
	y[idx+5] = rb.q[2];
	y[idx+6] = rb.q[3];

	y[idx+7] = rb.P[0];
	y[idx+8] = rb.P[1];
	y[idx+9] = rb.P[2];

	y[idx+10] = rb.L[0];
	y[idx+11] = rb.L[1];
	y[idx+12] = rb.L[2];
}

Simulation.prototype.normalize_quaterion = function(q){
  
	var s = q[0];
	var vx = q[1];
	var vy = q[2];
	var vz = q[3];

	var length = math.sqrt(s*s+vx*vx+vy*vy+vz*vz);
	return([s/length,vx/length,vy/length,vz/length]);
}

Simulation.prototype.array_to_state = function(rb, y, idx){

	rb.x[0] = y[idx+0];
	rb.x[1] = y[idx+1];
	rb.x[2] = y[idx+2];

	rb.q[0] = y[idx+3];
	rb.q[1] = y[idx+4];
	rb.q[2] = y[idx+5];
	rb.q[3] = y[idx+6];

	rb.P[0] = y[idx+7];
	rb.P[1] = y[idx+8];
	rb.P[2] = y[idx+9];
	
	rb.L[0] = y[idx+10];
	rb.L[1] = y[idx+11];
	rb.L[2] = y[idx+12];

	rb.v[0] = rb.P[0]/rb.mass;
	rb.v[1] = rb.P[1]/rb.mass;
	rb.v[2] = rb.P[2]/rb.mass;

	rb.R = this.quaterion_to_matrix(this.normalize_quaterion(rb.q));
	
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

Simulation.prototype.quaternion_product = function(a,b){
	//[s1,v1] * [s2,v2] = s1s2 - v1.v2, s1v2+s2v1 + v1xv2]
	var s1 = a[0];
	var v1 = [a[1],a[2],a[3]];
	var s2 = b[0];
	var v2 = [b[1],b[2],b[3]];

	var rs = s1*s2-this.dot_product(v1,v2);
	var rv = [0,0,0];
	
	var temp = this.cross_product(v1,v2);
	
	rv[0] = s1*v2[0] + s2*v1[0] + temp[0];	
	rv[1] = s1*v2[1] + s2*v1[1] + temp[1];
	rv[2] = s1*v2[2] + s2*v1[2] + temp[2];
	
	return [rs,rv[0],rv[1],rv[2]];
}

Simulation.prototype.dot_product = function(a,b){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

Simulation.prototype.cross_product = function(a,b){

  var a1 = a[0];
  var a2 = a[1];
  var a3 = a[2];
  var b1 = b[0];
  var b2 = b[1];
  var b3 = b[2];
  var first = a2*b3-a3*b2;
  var second = a3*b1-a1*b3;
  var third = a1*b2-a2*b1;
  return [first,second,third];

}

Simulation.prototype.compute_force_and_torque = function(t, rb){

  var f = [1,1,0];
  var r = new Array(3);
  r[0] = 0.002;
  r[1] = 0.003;
  r[2] = 0;
  
  rb.torque = this.cross_product(r,f);   
  
}

Simulation.prototype.DdtStateToArray = function(rb, xdot, idx){

	xdot[idx+0] = rb.v[0];
	xdot[idx+1] = rb.v[1];
	xdot[idx+2] = rb.v[2];
	
	var temp = [0,math.subset(rb.omega, math.index(0)), math.subset(rb.omega, math.index(1)), math.subset(rb.omega, math.index(2)) ];
	var qdot = this.quaternion_product(temp, rb.q);
	qdot[0] = 0.5*qdot[0];
	qdot[1] = 0.5*qdot[1];
	qdot[2] = 0.5*qdot[2];
	qdot[3] = 0.5*qdot[3];

	xdot[idx+3] = qdot[0];
	xdot[idx+4] = qdot[1];
	xdot[idx+5] = qdot[2];
	xdot[idx+6] = qdot[3];

	xdot[idx+7] = rb.force[0];
	xdot[idx+8] = rb.force[1];
	xdot[idx+9] = rb.force[2];

	xdot[idx+10] = rb.torque[0];
	xdot[idx+11] = rb.torque[1];
	xdot[idx+12] = rb.torque[2];
	
}

Simulation.prototype.Dxdt = function(t, x, xdot){

	this.array_to_bodies(x);
	for (var i=0; i<this.n_bodies; i++){
		this.compute_force_and_torque(t, this.rigid_bodies[i]);
		this.DdtStateToArray(this.rigid_bodies[i], xdot, i*this.state_size);
	}

}

Simulation.prototype.ode = function(x0, xFinal){

  xdot =  new Array(this.n_bodies * this.state_size);
  this.Dxdt(0, x0, xdot);
  for (var i=0;i<this.n_bodies * this.state_size; i++){
	xFinal[i] = x0[i] + xdot[i];	
  }
}

Simulation.prototype.draw = function(callback){
  var l1 = document.getElementById("l1");
  var l2 = document.getElementById("l2");
  var l3 = document.getElementById("l3");
  var l4 = document.getElementById("l4");

  var x1 = this.rigid_bodies[0].x[0];
  var y1 = this.rigid_bodies[0].x[1];
  
var w = 30;
  var x_1 = -w/2;
  var y_1 = -w/2;

  var x_2 = w/2;
  var y_2 = -w/2;

  var x_3 = w/2;
  var y_3 = w/2;

  var x_4 = -w/2;
  var y_4 = w/2;


  temp = math.matrix([[x_1,y_1,0]]);
  res = math.multiply(temp, this.rigid_bodies[0].R);

  l1.setAttribute("x1",  math.subset(res, math.index(0,0))+x1);
  l1.setAttribute("y1",  math.subset(res, math.index(0,1))+y1);

  l4.setAttribute("x2",  math.subset(res, math.index(0,0))+x1);
  l4.setAttribute("y2",  math.subset(res, math.index(0,1))+y1);

  temp = math.matrix([[x_2,y_2,0]]);
  res = math.multiply(temp, this.rigid_bodies[0].R);

  l1.setAttribute("x2",  math.subset(res, math.index(0,0))+x1);
  l1.setAttribute("y2",  math.subset(res, math.index(0,1))+y1);
 
  l2.setAttribute("x1",  math.subset(res, math.index(0,0))+x1);
  l2.setAttribute("y1",  math.subset(res, math.index(0,1))+y1);

  temp = math.matrix([[x_3,y_3,0]]);
  res = math.multiply(temp, this.rigid_bodies[0].R);
 
  l2.setAttribute("x2",  math.subset(res, math.index(0,0))+x1);
  l2.setAttribute("y2",  math.subset(res, math.index(0,1))+y1);
 
  l3.setAttribute("x1",  math.subset(res, math.index(0,0))+x1);
  l3.setAttribute("y1",  math.subset(res, math.index(0,1))+y1);

  temp = math.matrix([[x_4,y_4,0]]);
  res = math.multiply(temp, this.rigid_bodies[0].R);
   
  l3.setAttribute("x2",  math.subset(res, math.index(0,0))+x1);
  l3.setAttribute("y2",  math.subset(res, math.index(0,1))+y1);
 
  l4.setAttribute("x1",  math.subset(res, math.index(0,0))+x1);
  l4.setAttribute("y1",  math.subset(res, math.index(0,1))+y1);
  
}

Simulation.prototype.make_step = function(){
	for (var i=0;i<this.state_size*this.n_bodies;i++){
		this.x0[i] = this.xFinal[i];
	}
	this.ode(this.x0, this.xFinal);
	this.array_to_bodies(this.xFinal);
	this.draw();
	setTimeout(function(simulation){ simulation.make_step()}, 10, this);
}

Simulation.prototype.run_simulation = function(){

  this.x0 = new Array(this.n_bodies * this.state_size);
  this.xFinal = new Array(this.n_bodies * this.state_size);

  this.init_states();
  this.bodies_to_array(this.xFinal);

  this.make_step();
  
}
