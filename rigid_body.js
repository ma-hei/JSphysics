
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
  this.vertices = math.matrix([[1,1,0],[-1,1,0],[-1,-1,0],[1,-1,0]]);  
  this.translated_vertices;
  this.floor = false;
  
}


