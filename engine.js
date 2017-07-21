
var RigidBody = function(){
	
	this.mass = 1;
	this.I_body = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.I_body_inv = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);	
	this.x = [0,0,0];
	this.q = [0,0,0,1];
	this.P = [0,0,0];
	this.L = [0,0,0];
	this.I_inv = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.R = math.matrix([[1,0,0],[0,1,0],[0,0,1]]);
	this.v = [1,0,0];
	this.omega = [0,0,0];
	this.force = [0,0,0];
	this.torque = [0,0,0];
}


var Simulation = function(){

	this.n_bodies = 1;
	this.state_size = 13;
	this.rigid_bodies = new Array(this.n_bodies);
	this.stepsize = 1/60;
	this.time_step = 1/60;
	this.allowed_error = 0.001;
	this.steps_taken = 0;
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
	return([s/length,-1*vx/length,-1*vy/length,-1*vz/length]);
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

	debugger;
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

  if (t<0.5){
  	var f = [4,-9,0];
  	var r = new Array(3);
  	r[0] = 2;
  	r[1] = 1;
  	r[2] = 0;

  	rb.torque = this.cross_product(r,f);   
	  rb.force = f;
  	return;
  } else {
  	var f=[4,4,0];
  	var r=[2,-1,0];

  	rb.torque = this.cross_product(r,f);   
	  rb.force = f;
  }
  
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

Simulation.prototype.adapt_stepsize = function(x0, t){
	temp = new Array(this.n_bodies * this.state_size);
	debugger;
	this.Dxdt(t, x0, temp);
	a = new Array(this.n_bodies * this.state_size);
	for (var i=0; i<this.n_bodies * this.state_size; i++){
		a[i] = x0[i]+this.stepsize*temp[i];
	}

	b = new Array(this.n_bodies * this.state_size);
	for (var i=0; i<this.n_bodies * this.state_size; i++){
		b[i] = x0[i]+0.5*this.stepsize*temp[i];
	}
	debugger;
	this.Dxdt(t+(this.stepsize), b, temp);
	for (var i=0;i<this.n_bodies * this.state_size; i++){
		b[i] = b[i]+0.5*this.stepsize*temp[i];
	}

	var error = 0;
	for (var i=0;i<this.n_bodies * this.state_size; i++){
		error = error + (a[i]-b[i])*(a[i]-b[i]);
	}
	error = Math.sqrt(error);
	//console.log("error: "+error);
	var new_stepsize = Math.sqrt(this.allowed_error/error) * this.stepsize;
	this.stepsize = Math.min(this.time_step,new_stepsize);
}

Simulation.prototype.euler_step = function (x0, xFinal, t, t_end, stepsize){
	xdot = new Array(this.n_bodies * this.state_size);
	this.Dxdt(t, x0, xdot);
	for (var i=0;i<this.n_bodies*this.state_size;i++){
		xFinal[i] = x0[i] + stepsize*xdot[i];
	}
}

Simulation.prototype.euler_step_2 = function(x0, xFinal, t, t_end){
	var current_time = t;
	for (var i=0;i<this.n_bodies*this.state_size;i++){
		xFinal[i] = x0[i];
	}
	while(current_time<t_end){
		this.adapt_stepsize(x0,t);
		var step = Math.min(t_end-current_time, this.stepsize);
		xdot = new Array(this.n_bodies * this.state_size);
		this.Dxdt(t, x0, xdot);
		for (var i=0;i<this.n_bodies*this.state_size;i++){
			xFinal[i] = xFinal[i] + step*xdot[i];
		}
		current_time = current_time + step;
	}
}

Simulation.prototype.runge_katta = function(x0, xFinal, current_time, stepsize){
    var array_size = this.n_bodies * this.state_size;
    var temp = new Array(array_size);
		var k1 = new Array(array_size);
    var k2 = new Array(array_size);
    var k3 = new Array(array_size);
    var k4 = new Array(array_size); 
  	this.Dxdt(current_time, x0, k1);
    for (var i=0; i<array_size; i++){
      k1[i] = stepsize*k1[i];
    }

    // compute k2
    for (var i=0;i<array_size;i++){
      temp[i] = x0[i] + k1[i]/2;
    }
    this.Dxdt(current_time+stepsize/2, temp, k2);
    for (var i=0; i<array_size; i++){
      k2[i] = stepsize*k2[i];
    }

    // compute k3
    for (var i=0;i<array_size;i++){
      temp[i] = x0[i] + k2[i]/2;
    }
    this.Dxdt(current_time+stepsize/2, temp, k3);
    for (var i=0; i<array_size; i++){
      k3[i] = stepsize*k3[i];
    }

    // compute k4
    for (var i=0;i<array_size;i++){
      temp[i] = x0[i] + k3[i];
    }
    this.Dxdt(current_time+stepsize, temp, k4);
    for (var i=0; i<array_size; i++){
      k4[i] = stepsize*k4[i];
    }

    // compute xFinal
    for (var i=0;i<array_size;i++){
      xFinal[i] = x0[i] + (1/6)*k1[i] + (1/3)*k2[i] + (1/3)*k3[i] + (1/6)*k4[i];
    }

}

Simulation.prototype.ode = function(x0, xFinal,t, t_end){

  this.runge_katta(x0, xFinal, t, this.time_step);
  this.compare_error(x0,t);

}

Simulation.prototype.compare_error = function(x0, t){

  var array_size = this.n_bodies * this.state_size;
  var temp = new Array(array_size);
  this.runge_katta_2(x0, temp, t, this.time_step);
  
  var temp2 = new Array(array_size);
  this.runge_katta_2(x0, temp2, t, this.time_step/2);
  var temp3 = new Array(array_size);
  this.runge_katta_2(temp2, temp3, t+this.time_step/2, this.time_step/2);

  var error = 0;
  for (var i=0; i<array_size; i++){
    error += (temp3[i] - temp[i])*(temp3[i]-temp[i]);
  }
  error = Math.sqrt(error);

  this.euler_step(x0, temp, t, 0, this.time_step);

  this.euler_step(x0, temp2, t, 0, this.time_step/2);
  this.euler_step(temp2, temp3, t+this.time_step/2, 0, this.time_step/2);

  var error2 = 0;
  for (var i=0; i<array_size; i++){
    error2 += (temp3[i] - temp[i])*(temp3[i]-temp[i]);
  }
  error2 = Math.sqrt(error2);

  console.log("error rk: "+error+", error es: "+error2);


}

Simulation.prototype.draw = function(){

  var x1 = this.rigid_bodies[0].x[0];
  var y1 = this.rigid_bodies[0].x[1];
	
	var temp = this.rigid_bodies[0].R;
	var temp = mat4.create();
	temp[0] = math.subset(this.rigid_bodies[0].R, math.index(0,0));
	temp[1] = math.subset(this.rigid_bodies[0].R, math.index(0,1));
	temp[2] = math.subset(this.rigid_bodies[0].R, math.index(0,2));
	temp[4] = math.subset(this.rigid_bodies[0].R, math.index(1,0));
	temp[5] = math.subset(this.rigid_bodies[0].R, math.index(1,1));
	temp[6] = math.subset(this.rigid_bodies[0].R, math.index(1,2));
	temp[11] = 1;
	temp[15] = 1;

	gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  mat4.perspective(45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0, pMatrix);

  mat4.identity(mvMatrix);

  mat4.translate(mvMatrix, [x1, y1, -80.0]);

	mvPushMatrix();
	
	mat4.multiply(mvMatrix, temp, mvMatrix);
  gl.bindBuffer(gl.ARRAY_BUFFER, squareVertexPositionBuffer);
  gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, squareVertexPositionBuffer.itemSize, gl.FLOAT, false, 0, 0);
  setMatrixUniforms();
  gl.drawArrays(gl.TRIANGLE_STRIP, 0, squareVertexPositionBuffer.numItems);
	mvPopMatrix();

}

Simulation.prototype.make_step = function(t){
	for (var i=0;i<this.state_size*this.n_bodies;i++){
		this.x0[i] = this.xFinal[i];
	}
	this.ode(this.x0, this.xFinal, t, t+simulation.time_step);
	this.array_to_bodies(this.xFinal);
	this.draw();
	this.steps_taken = this.steps_taken+1;
	setTimeout(function(simulation,t){ simulation.make_step(t)}, simulation.time_step*1000, this,t+simulation.time_step);
}

Simulation.prototype.run_simulation = function(){
	
	webGLStart();

  this.x0 = new Array(this.n_bodies * this.state_size);
  this.xFinal = new Array(this.n_bodies * this.state_size);

  this.init_states();
  this.bodies_to_array(this.xFinal);

  this.make_step(0);
  
}
