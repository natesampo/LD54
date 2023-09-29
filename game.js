class Vector {
	constructor(x, y) {
		this.x = x;
		this.y = y;
	}

	translate(vector) {
		this.x += vector.x;
		this.y += vector.y;
	}

	addVector(vector) {
		this.x += vector.x;
		this.y += vector.y;
	}

	subtractVector(vector) {
		this.x -= vector.x;
		this.y -= vector.y;
	}

	addScalar(scalar) {
		this.x += scalar;
		this.y += scalar;
	}

	multiplyScalar(scalar) {
		this.x *= scalar;
		this.y *= scalar;
	}

	divideScalar(scalar) {
		this.x /= scalar;
		this.y /= scalar;
	}

	differenceVector(vector) {
		return new Vector(this.x - vector.x, this.y - vector.y);
	}

	dotProduct(vector) {
		return this.x * vector.x + this.y * vector.y;
	}

	vectorCrossProduct(vector) {
		return this.x * vector.y - this.y * vector.x;
	}

	scalarCrossProduct(scalar) {
		return new Vector(scalar * this.y, -scalar * this.x);
	}

	magnitude() {
		return Math.sqrt(this.x * this.x + this.y * this.y);
	}

	rotateAboutPoint(vector, angle) {
		let cos = Math.cos(angle);
		let sin = Math.sin(angle);

		let newX = (cos * (this.x - vector.x)) - (sin * (this.y - vector.y)) + vector.x;
		let newY = (cos * (this.y - vector.y)) + (sin * (this.x - vector.x)) + vector.y;

		this.x = newX;
		this.y = newY;
	}

	negate() {
		this.x = -this.x;
		this.y = -this.y;
	}

	normalize() {
		let distance = Math.sqrt(this.x * this.x + this.y * this.y);
		if (distance > 0) {
			this.x /= distance;
			this.y /= distance;
		}
	}

	set(x, y) {
		this.x = x;
		this.y = y;
	}

	copy() {
		return new Vector(this.x, this.y);
	}
}
let gravity = new Vector(0, 4);

class Link {
	constructor(position) {
		this.position = position;

		this.size = 20;
		this.angle = 0;
		this.red = 0;
		this.green = 0;
		this.blue = 0;
		this.alpha = 1;

		this.body1 = null;
		this.body2 = null;
	}

	findLinks(bodies) {
		for (var i=bodies.length-1; i>=0; i--) {
			if (bodies[i].isInside(this.position)) {
				if (this.body1 == null) {
					this.body1 = bodies[i];
				} else {
					this.body2 = bodies[i];
					break;
				}
			}
		}
	}

	getPartner(body) {
		if (body == this.body1) {
			return this.body2;
		}
		if (body == this.body2) {
			return this.body1;
		}

		return null;
	}

	translate(vector) {
		this.position.translate(vector);
	}

	setXY(vector) {
		this.position.set(vector.x, vector.y);
	}

	rotate(rotation) {
		this.angle += rotation;
		this.angle = this.angle % (2*Math.PI);
		if (this.angle < 0) {
			this.angle += 2*Math.PI;
		}
	}

	isInside(vector) {
		return getDistance(this.position, vector) < this.size;
	}

	render(context, camera) {
	}

	getName() {
		return 'Link';
	}
}

class Rivet extends Link {
	constructor(position) {
		super(position);
		this.blue = 255;
		this.compositeObject = null;

		this.relativeAngle = null;
		this.relativePosition = null;
		this.relativeTheta = null;
		this.initialTheta = null;
	}

	findLinks(bodies) {
		super.findLinks(bodies);

		if (this.body1 != null) {
			this.body1.rivets.push(this);

			if (this.body1.compositeObject == null) {
				if (this.body2 == null || this.body2.compositeObject == null) {
					new CompositeObject([this.body1, this.body2], [this]);
				} else if (this.body2 != null && this.body2.compositeObject != null) {
					this.body2.compositeObject.addObject(this.body1);
					this.body2.compositeObject.addLink(this);
				}
			} else {
				if (this.body2 == null || this.body2.compositeObject == null) {
					this.body1.compositeObject.addObject(this.body2);
					this.body1.compositeObject.addLink(this);
				} else if (this.body2 != null && this.body2.compositeObject != null) {
					this.body1.compositeObject.addLink(this);
					this.body1.compositeObject.merge(this.body2.compositeObject);
				}
			}
		}
		if (this.body2 != null) {
			this.body2.rivets.push(this);
		}
	}

	updatePosition() {
		if (this.body1 != null) {
			if (this.relativeAngle == null) {
				this.relativeAngle = this.body1.angle - this.angle;
				this.relativePosition = getDistance(this.body1.vertices[0], this.position);
				this.relativeTheta = getAngle(this.body1.vertices[0], this.position);
				this.initialTheta = this.body1.angle;
			} else {
				this.rotate((this.body1.angle - this.relativeAngle) - this.angle);
				this.setXY(new Vector(this.body1.vertices[0].x + this.relativePosition * Math.cos(this.body1.angle - this.initialTheta + this.relativeTheta), this.body1.vertices[0].y + this.relativePosition * Math.sin(this.body1.angle - this.initialTheta + this.relativeTheta)));
			}
		}
	}

	delete() {
		if (this.body1 != null) {
			remove(this.body1.rivets, this);
		}
		if (this.body2 != null) {
			remove(this.body2.rivets, this);
		}

		this.compositeObject.removeLink(this);

		this.body1 = null;
		this.body2 = null;
	}


	render(context, camera) {
		let canvas = context.canvas;

		let newCanvas = document.createElement('canvas');
		newCanvas.width = this.size;
		newCanvas.height = this.size;
		let newContext = newCanvas.getContext('2d');

		newContext.strokeStyle = 'rgba(0, 0, 0, ' + ((this.body1 == null && this.body2 == null) ? 0.5 : 1) + ')';
		newContext.lineWidth = 5;
		newContext.beginPath();
		newContext.moveTo(this.size/2 - (this.size/2)*Math.cos(this.angle + Math.PI/4), this.size/2 - (this.size/2)*Math.sin(this.angle + Math.PI/4));
		newContext.lineTo(this.size/2 + (this.size/2)*Math.cos(this.angle + Math.PI/4), this.size/2 + (this.size/2)*Math.sin(this.angle + Math.PI/4));
		newContext.moveTo(this.size/2 + (this.size/2)*Math.sin(this.angle + Math.PI/4), this.size/2 - (this.size/2)*Math.cos(this.angle + Math.PI/4));
		newContext.lineTo(this.size/2 - (this.size/2)*Math.sin(this.angle + Math.PI/4), this.size/2 + (this.size/2)*Math.cos(this.angle + Math.PI/4));
		newContext.stroke();
		newContext.closePath();
		newContext.beginPath();
		newContext.moveTo(this.size/2 - (this.size/2 - 1)*Math.cos(this.angle + Math.PI/4), this.size/2 - (this.size/2 - 1)*Math.sin(this.angle + Math.PI/4));
		newContext.lineTo(this.size/2 + (this.size/2 - 1)*Math.cos(this.angle + Math.PI/4), this.size/2 + (this.size/2 - 1)*Math.sin(this.angle + Math.PI/4));
		newContext.moveTo(this.size/2 + (this.size/2 - 1)*Math.sin(this.angle + Math.PI/4), this.size/2 - (this.size/2 - 1)*Math.cos(this.angle + Math.PI/4));
		newContext.lineTo(this.size/2 - (this.size/2 - 1)*Math.sin(this.angle + Math.PI/4), this.size/2 + (this.size/2 - 1)*Math.cos(this.angle + Math.PI/4));
		newContext.lineWidth = 3;
		newContext.globalCompositeOperation = 'destination-out';
		newContext.stroke();
		newContext.globalCompositeOperation = 'source-over';
		newContext.strokeStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + this.alpha + ')';
		newContext.stroke();
		newContext.closePath();

		context.drawImage(newCanvas, this.position.x - this.size/2 - camera.x, this.position.y - this.size/2 - camera.y);
	}

	getName() {
		return 'Rivet';
	}
}

class Hinge extends Link {
	constructor(position) {
		super(position);
		this.strength = 1000;
		this.motorSpeed = 0;

		this.body1Distance = null;
		this.body2Distance = null;
		this.body1Theta = null;
		this.body2Theta = null;
		this.body1InitialAngle = null;
		this.body2InitialAngle = null;
	}

	render(context, camera) {
		let canvas = context.canvas;

		context.strokeStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + ((this.body1 == null && this.body2 == null) ? 0.5 : 1) + ')';
		context.lineWidth = 3;
		context.beginPath();
		context.arc(this.position.x - camera.x, this.position.y - camera.y, this.size/2, 0, 2*Math.PI, false);
		context.stroke();
		context.closePath();
	}

	initPoints() {
		let body1COM = this.body1.findCenterOfMass();
		let body2COM = this.body2.findCenterOfMass();

		this.body1Distance = getDistance(body1COM, this.position);
		this.body2Distance = getDistance(body2COM, this.position);
		this.body1Theta = getAngle(body1COM, this.position);
		this.body2Theta = getAngle(body2COM, this.position);
		this.body1InitialAngle = this.body1.angle;
		this.body2InitialAngle = this.body2.angle;
	}

	applyImpulse() {
		let body1COM = this.body1.findCenterOfMass();
		let body2COM = this.body2.findCenterOfMass();

		let body1AnchorPoint = new Vector(body1COM.x + this.body1Distance * Math.cos(this.body1.angle - this.body1InitialAngle + this.body1Theta), body1COM.y + this.body1Distance * Math.sin(this.body1.angle - this.body1InitialAngle + this.body1Theta));
		let body2AnchorPoint = new Vector(body2COM.x + this.body2Distance * Math.cos(this.body2.angle - this.body2InitialAngle + this.body2Theta), body2COM.y + this.body2Distance * Math.sin(this.body2.angle - this.body2InitialAngle + this.body2Theta));

		this.position.set((body1AnchorPoint.x + body2AnchorPoint.x)/2, (body1AnchorPoint.y + body2AnchorPoint.y)/2);

		let directionVector = body1AnchorPoint.differenceVector(body2AnchorPoint);
		let separation = directionVector.magnitude();
		directionVector.multiplyScalar(this.strength);

		body1AnchorPoint.subtractVector(body1COM);
		body2AnchorPoint.subtractVector(body2COM);

		this.body2.applyImpulse(directionVector, body2AnchorPoint);
		directionVector.negate();
		this.body1.applyImpulse(directionVector, body1AnchorPoint);
		
		this.spinMeRightRound();
	}

	spinMeRightRound() {
		if (this.body1 != null && this.body2 != null) {
			if (this.body1.getMass() == 0) {
				this.body2.torque -= (this.motorSpeed/100);
			} else if (this.body2.getMass() == 0) {
				this.body1.torque += (this.motorSpeed/100);
			} else {
				this.body1.torque += (this.motorSpeed/100) * (1-(this.body1.getMass()/(this.body2.getMass() + this.body1.getMass())));
				this.body2.torque -= (this.motorSpeed/100) * (1-(this.body2.getMass()/(this.body2.getMass() + this.body1.getMass())));
			}
		} else {
			if (this.body1 != null) {
				this.body1.torque += (this.motorSpeed/100);
			} else if (this.body2 != null) {
				this.body2.torque -= (this.motorSpeed/100);
			}
		}
	}

	findLinks(bodies) {
		super.findLinks(bodies);

		if (this.body1 != null) {
			this.body1.hinges.push(this);
		}
		if (this.body2 != null) {
			this.body2.hinges.push(this);
		}
	}

	delete() {
		if (this.body1 != null) {
			remove(this.body1.hinges, this);
		}
		if (this.body2 != null) {
			remove(this.body2.hinges, this);
		}

		this.body1 = null;
		this.body2 = null;
	}

	getName() {
		return 'Hinge';
	}
}

class Body {
	constructor(vertices) {
		this.vertices = vertices;
		this.id = 0;
		this.normals = [];
		this.density = 1;
		this.area = 1;
		this.angle = 0;
		this.mass = 1;
		this.invMass = 1;
		this.momentOfInertia = 1;
		this.invMomentOfInertia = 1/this.momentOfInertia;
		this.restitution = 0.2;
		this.friction = 0.1;
		this.red = 255;
		this.green = 0;
		this.blue = 0;
		this.alpha = 1;
		this.rivetedToBackground = false;
		this.rivets = [];
		this.compositeObject = null;
		this.hinges = [];

		this.velocity = new Vector(0, 0);
		this.force = new Vector(0, 0);

		this.angularVelocity = 0;
		this.torque = 0;

		this.computeArea();
		this.computeMass();

		if (this.vertices.length > 1) {
			this.computeNormals();
		}
	}

	applyForce(vector) {
		this.force.addVector(vector);
	}

	applyImpulse(impulse, contactVector) {
		if (this.compositeObject == null) {
			let tempImpulse = impulse.copy();
			tempImpulse.multiplyScalar(this.getInvMass());
			this.velocity.addVector(tempImpulse);

			let deltaAngularVelocity = this.getInvMomentOfInertia() * contactVector.vectorCrossProduct(impulse);
			this.angularVelocity += deltaAngularVelocity;
		} else {
			this.compositeObject.applyImpulse(impulse, contactVector);
		}
	}

	translateComposite(vector) {
		if (this.compositeObject == null) {
			this.translate(vector);
		} else {
			this.compositeObject.translate(vector);
		}
	}

	translate(vector) {
		for(var i=0; i<this.vertices.length; i++) {
			this.vertices[i].translate(vector);
		}
	}

	rotate(rotation) {
		if (this.vertices.length > 0) {
			if (this.vertices.length > 1 || this.compositeObject != null) {
				if (this.hinges.length > 0 && this.hinges[0].getPartner(this) == null) {
					for (var i=0; i<this.vertices.length; i++) {
						this.vertices[i].rotateAboutPoint(this.hinges[0].position, rotation);
					}
				} else {
					let centerOfMass = this.findCenterOfMass();
					for (var i=0; i<this.vertices.length; i++) {
						this.vertices[i].rotateAboutPoint(centerOfMass, rotation);
					}
				}
			}

			this.angle += rotation;
			this.angle = this.angle % (2*Math.PI);
			if (this.angle < 0) {
				this.angle += 2*Math.PI;
			}

			if (this.vertices.length > 1) {
				this.computeNormals();
			}
		}
	}

	applyHingeImpulses() {
		for (var i=0; i<this.hinges.length; i++) {
			if (this.hinges[i].body1Distance == null) {
				if (this.hinges[i].getPartner(this) == null) {
					this.setStatic();
					if (this.compositeObject != null) {
						this.compositeObject.calculateMass();
					}
					this.hinges[i].spinMeRightRound();
				} else {
					this.hinges[i].initPoints();
				}
			}

			if (this.getMass() > 0) {
				this.hinges[i].applyImpulse();
			}
		}
	}

	getFarthestVertex(vector) {
		let farthestProjection = null;
		let farthestVertex = null;

		for (var i=0; i<this.vertices.length; i++) {
			let projection = this.vertices[i].dotProduct(vector);

			if (farthestProjection == null || projection > farthestProjection) {
				farthestProjection = projection;
				farthestVertex = this.vertices[i];
			}
		}

		return farthestVertex;
	}

	getLeastPenetration(body) {
		let mostDistance = null;
		let bestFace = null;

		for (var i=0; i<this.vertices.length; i++) {
			let normal = this.normals[i].copy();
			normal.negate();
			let farthestVertex = body.getFarthestVertex(normal).copy();
			farthestVertex.subtractVector(this.vertices[i]);
			normal.negate();

			let distance = normal.dotProduct(farthestVertex);
			if (mostDistance == null || distance > mostDistance) {
				mostDistance = distance;
				bestFace = i;
			}
		}

		return [mostDistance, bestFace];
	}

	findIncidentFace(incidentBody, referenceFace) {
		let incidentFace = null;
		let minimumDot = null;
		for (var i=0; i<incidentBody.vertices.length; i++) {
			let dot = this.normals[referenceFace].dotProduct(incidentBody.normals[i]);
			if (minimumDot == null || dot < minimumDot) {
				minimumDot = dot;
				incidentFace = i;
			}
		}

		return incidentFace;
	}

	getMass() {
		if (this.compositeObject == null) {
			return this.mass;
		}
		
		return this.compositeObject.getMass();
	}

	getInvMass() {
		if (this.compositeObject == null) {
			return this.invMass;
		}
		
		return this.compositeObject.getInvMass();
	}

	getInvMomentOfInertia() {
		if (this.compositeObject == null) {
			return this.invMomentOfInertia;
		}
		
		return this.compositeObject.getInvMomentOfInertia();
	}

	findCenterOfMass() {
		if (this.compositeObject == null) {
			return this.findCenterOfMassForComposite();
		}
		
		return this.compositeObject.findCenterOfMass();
	}

	findCenterOfMassForComposite() {
		let centerOfMass = new Vector(0, 0);
		for (var i=0; i<this.vertices.length; i++) {
			centerOfMass.addVector(this.vertices[i]);
		}
		centerOfMass.divideScalar(this.vertices.length);

		return centerOfMass;
	}

	computeArea() {
		this.area = 1;
	}

	computeMass() {
		this.mass = this.density * this.area;
		this.invMass = (this.mass == 0) ? 0 : 1/this.mass;

		if (this.compositeObject != null) {
			this.compositeObject.calculateMass();
		}
	}

	computeNormals() {
		this.normals = [];
		for (var i=0; i<this.vertices.length; i++) {
			let normalVector = this.vertices[(i+1) % this.vertices.length].copy();
			normalVector.subtractVector(this.vertices[i]);
			normalVector.normalize();

			let oldX = normalVector.x;
			normalVector.x = normalVector.y;
			normalVector.y = -oldX;

			this.normals.push(normalVector);
		}
	}

	setStatic() {
		this.mass = 0;
		this.invMass = 0;
		this.momentOfInertia = 0;
		this.invMomentOfInertia = 0;
	}

	isInside(vector) {
		return false;
	}

	delete() {
		for (var i=this.rivets.length-1; i>=0; i--) {
			this.rivets[i].delete();
		}

		this.rivets = [];
	}

	tick(dt) {
		if (this.getMass() > 0 && !this.rivetedToBackground) {
			let dtForce = this.force.copy();
			dtForce.multiplyScalar(this.getInvMass());
			dtForce.addVector(gravity);
			dtForce.multiplyScalar(dt);
			this.velocity.addVector(dtForce);
			this.translate(this.velocity);

			if (Math.abs(this.angularVelocity) < 0.1) {
				this.angularVelocity += this.torque * dt;
			}
			if (this.angularVelocity != 0) {
				this.rotate(this.angularVelocity);
			}

			let airResistance = 0.9999;
			this.velocity.multiplyScalar(airResistance);
			this.angularVelocity *= airResistance;

			for (var i=0; i<this.rivets.length; i++) {
				if (this.rivets[i].body1 == null || this.rivets[i].body1 == this) {
					this.rivets[i].updatePosition();
				}
			}
		} else if (this.getMass() == 0) {
			this.angularVelocity += this.torque * dt;
			if (this.angularVelocity != 0) {
				this.rotate(this.angularVelocity);
			}
		}

		this.force.set(0, 0);
		this.torque = 0;
	}

	render(context, camera) {
		let canvas = context.canvas;
		context.strokeStyle = 'rgba(0, 0, 0, 1)';
		context.fillStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + this.alpha + ')';
		context.lineWidth = 3;
		context.beginPath();
		context.moveTo(this.vertices[0].x - camera.x, this.vertices[0].y - camera.y);
		for (var i=1; i<this.vertices.length; i++) {
			context.lineTo(this.vertices[i].x - camera.x, this.vertices[i].y - camera.y);
		}
		context.lineTo(this.vertices[0].x - camera.x, this.vertices[0].y - camera.y);
		context.fill();
		context.stroke();
		context.closePath();
	}

	getName() {
		return 'Body';
	}
}

class Circle extends Body {
	constructor(center, radius) {
		super([center]);
		this.radius = radius;
		this.computeArea();
		this.computeMass();
		this.momentOfInertia = 4000000;
		this.invMomentOfInertia = 1/this.momentOfInertia;
	}

	computeArea() {
		if (this.radius) {
			this.area = Math.PI * this.radius * this.radius;
		}
	}

	isInside(vector) {
		return getDistance(vector, this.vertices[0]) < this.radius;
	}

	render(context, camera) {
		let canvas = context.canvas;

		context.strokeStyle = 'rgba(0, 0, 0, 1)';
		context.fillStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + this.alpha + ')';
		context.lineWidth = 3;
		context.beginPath();
		context.arc(this.vertices[0].x - camera.x, this.vertices[0].y - camera.y, this.radius, 0, 2*Math.PI, false);
		context.fill();
		context.moveTo(this.vertices[0].x - camera.x, this.vertices[0].y - camera.y);
		context.lineTo(this.vertices[0].x - camera.x + Math.cos(this.angle) * this.radius, this.vertices[0].y - camera.y + Math.sin(this.angle) * this.radius);
		context.stroke();
		context.closePath();
	}

	getName() {
		return 'Circle';
	}
}

class Rectangle extends Body {
	constructor(vertices) {
		super(vertices);
		this.momentOfInertia = 100000000;
		this.invMomentOfInertia = 1/this.momentOfInertia;
	}

	computeArea() {
		this.area = getDistance(this.vertices[0], this.vertices[1]) * getDistance(this.vertices[1], this.vertices[2]);
	}

	isInside(vector) {
		let tempArea = areaOfTriangle(this.vertices[0], vector, this.vertices[3]) + areaOfTriangle(this.vertices[3], vector, this.vertices[2]) + 
						areaOfTriangle(this.vertices[2], vector, this.vertices[1]) + areaOfTriangle(vector, this.vertices[1], this.vertices[0]);
		return Math.floor(tempArea) <= Math.floor(this.area);
	}

	getName() {
		return 'Rectangle';
	}
}

class Camera {
	constructor() {
		this.x = 20;
		this.y = 20;
		this.zoom = 1;

		this.xVel = 0;
		this.yVel = 0;
	}
}

class Game {
	constructor() {
		this.bodies = [];
		this.currid = 0;

		this.levels = [[new Rectangle([new Vector(200, 1000), new Vector(1600, 1000), new Vector(1600, 1100), new Vector(200, 1100)])]];
		this.levels[0][0].setStatic();

		this.level = 0;

		this.ticksPerSecond = 120;

		this.mousePosition = new Vector(0, 0);
		this.oldMousePosition = new Vector(0, 0);
		this.dragging = null;
		this.placing = null;
		this.placingDrag = false;
		this.draggingSlider = null;
		this.startTime = null;
		this.tickID = 0;
		this.dt = 1/this.ticksPerSecond;
		this.collisions = [];
		this.contextMenu = null;

		this.buildMenuY = 200;
		this.buildMenuWidth = 50;
		this.buildMenuHeight = 400;
		this.buildMenuCurveSize = 24;
		this.buildMenuIconSize = 48;

		this.contextMenuWidth = 300;
		this.contextMenuHeight = 400;
		this.contextMenuCurveSize = 16;
		this.contextMenuSliderWidth = 0.75 * this.contextMenuWidth;
		this.contextMenuSliderThickness = 3;
		this.contextMenuSliderGrabberWidth = 10;
		this.contextMenuSliderGrabberHeight = 14;
		this.contextMenuSliderGap = 45;
		this.bodySliders = ['angle', 'density', 'restitution', 'friction', 'red', 'green', 'blue', 'alpha'];
		this.bodySlidersMin = [0, 0.1, 0, 0, 0, 0, 0, 0];
		this.bodySlidersMax = [Math.round(2*Math.PI*100)/100, 10, 1, 1, 255, 255, 255, 1];
		this.rivetSliders = ['size', 'angle', 'red', 'green', 'blue', 'alpha'];
		this.rivetSlidersMin = [5, 0, 0, 0, 0, 0];
		this.rivetSlidersMax = [100, Math.round(2*Math.PI*100)/100, 255, 255, 255, 1];
		this.hingeSliders = ['motorSpeed', 'size', 'red', 'green', 'blue'];
		this.hingeSlidersMin = [-11, 5, 0, 0, 0];
		this.hingeSlidersMax = [11, 100, 255, 255, 255];

		this.canvas = document.createElement('canvas');
		this.canvas.width = window.innerWidth;
		this.canvas.height = window.innerHeight;
		this.canvas.oncontextmenu = function() {return false;};
		this.context = this.canvas.getContext('2d');
		this.context.imageSmoothingEnabled = false;
		this.context.mozImageSmoothingEnabled = false;
		this.context.webkitImageSmoothingEnabled = false;

		this.inputs = [];
		this.camera = new Camera();

		this.setupListeners();

		this.pause = true;
		this.canInteract = true;
		this.debug = false;

		this.loadLevel(this.levels[this.level]);
	}

	loadLevel(level) {
		this.contextMenu = null;
		this.dragging = null;
		this.placing = null;
		this.placingDrag = false;

		this.bodies = [];
		for (var i=0; i<level.length; i++) {
			this.addBody(level[i]);
		}

		this.camera.x = 0;
		this.camera.y = 0;
	}

	setupListeners() {
		let gameForListeners = this;
		addMouseDownListener(function(which, eventX, eventY) {
			let eventPosition = new Vector(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y);
			switch(which) {
				case 1:
					let inContextMenu = false;
					if (gameForListeners.contextMenu != null) {
						let contextX = gameForListeners.contextMenu[1].x;
						let contextY = gameForListeners.contextMenu[1].y;
						let width = gameForListeners.contextMenuWidth;
						let height = gameForListeners.contextMenuHeight;
						let curveSize = gameForListeners.contextMenuCurveSize;

						if (eventX < contextX || eventX > contextX + 2*curveSize + width || eventY < contextY - 2*curveSize - height || eventY > contextY) {
							gameForListeners.contextMenu = null;
						} else {
							inContextMenu = true;
							let grabbedSlider = false;
							if (eventX >= gameForListeners.contextMenu[1].x + gameForListeners.contextMenuCurveSize + (gameForListeners.contextMenuWidth - gameForListeners.contextMenuSliderWidth)/2 &&
								eventX <= gameForListeners.contextMenu[1].x + gameForListeners.contextMenuCurveSize + (gameForListeners.contextMenuWidth - gameForListeners.contextMenuSliderWidth)/2 + gameForListeners.contextMenuSliderWidth) {
								let sliders = [];
								if (gameForListeners.contextMenu[0] instanceof Body) {
									sliders = gameForListeners.bodySliders;
								} else if (gameForListeners.contextMenu[0] instanceof Rivet) {
									sliders = gameForListeners.rivetSliders;
								} else if (gameForListeners.contextMenu[0] instanceof Hinge) {
									sliders = gameForListeners.hingeSliders;
								}

								for (var i=0; i<sliders.length; i++) {
									let sliderHeight = gameForListeners.contextMenuHeight - gameForListeners.contextMenuCurveSize - 22 - gameForListeners.contextMenuSliderGap * i;
									if (eventY >= gameForListeners.contextMenu[1].y - sliderHeight - gameForListeners.contextMenuSliderGrabberHeight - gameForListeners.contextMenuSliderThickness &&
										eventY <= gameForListeners.contextMenu[1].y - sliderHeight + gameForListeners.contextMenuSliderGrabberHeight + gameForListeners.contextMenuSliderThickness) {
										gameForListeners.draggingSlider = i;
										grabbedSlider = true;
										break;
									}
								}
							}

							if (!grabbedSlider) {
								let top = gameForListeners.contextMenu[1].y - 2*gameForListeners.contextMenuCurveSize - gameForListeners.contextMenuHeight
								if (eventX < gameForListeners.contextMenu[1].x + 2*gameForListeners.contextMenuCurveSize + gameForListeners.contextMenuWidth - 12 &&
									eventX > gameForListeners.contextMenu[1].x + 2*gameForListeners.contextMenuCurveSize + gameForListeners.contextMenuWidth - 12 - gameForListeners.contextMenuCurveSize && 
									eventY < top + 8 + gameForListeners.contextMenuCurveSize && 
									eventY > top + 8) {
									if (confirm('Delete this object?')) {
										if (gameForListeners.contextMenu[0] instanceof Body) {
											gameForListeners.removeBody(gameForListeners.contextMenu[0]);
										} else if (gameForListeners.contextMenu[0] instanceof Rivet) {
											gameForListeners.contextMenu[0].delete();
											gameForListeners.contextMenu = null;
										} else if (gameForListeners.contextMenu[0] instanceof Hinge) {
											gameForListeners.contextMenu[0].delete();
											gameForListeners.contextMenu = null;
										}
									}
								}
							}
						}
					}

					if (gameForListeners.canInteract && !inContextMenu) {
						if (gameForListeners.placing != null) {
							gameForListeners.placingDrag = true;
						} else {
							if (eventX < gameForListeners.buildMenuWidth + gameForListeners.buildMenuCurveSize && 
								eventY < gameForListeners.buildMenuY + 2*gameForListeners.buildMenuCurveSize + gameForListeners.buildMenuHeight && eventY > gameForListeners.buildMenuY) {
								if (eventX < (gameForListeners.buildMenuWidth + gameForListeners.buildMenuCurveSize/1.25 + gameForListeners.buildMenuIconSize)/2 && 
									eventX > (gameForListeners.buildMenuWidth + gameForListeners.buildMenuCurveSize/1.25 - gameForListeners.buildMenuIconSize)/2) {
									for (var i=0; i<4; i++) {
										if (eventY < gameForListeners.buildMenuY + gameForListeners.buildMenuCurveSize + i*(gameForListeners.buildMenuIconSize + 1.5*gameForListeners.buildMenuCurveSize) + gameForListeners.buildMenuIconSize && 
											eventY > gameForListeners.buildMenuY + gameForListeners.buildMenuCurveSize + i*(gameForListeners.buildMenuIconSize + 1.5*gameForListeners.buildMenuCurveSize)) {
											switch(i) {
												case 0:
													gameForListeners.placing = new Circle(eventPosition.copy(), 50);
													gameForListeners.placing.alpha = 0.5;
													break;
												case 1:
													gameForListeners.placing = new Rectangle([new Vector(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y),
																				new Vector(eventX + gameForListeners.camera.x + 200, eventY + gameForListeners.camera.y),
																				new Vector(eventX + gameForListeners.camera.x + 200, eventY + gameForListeners.camera.y + 100),
																				new Vector(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y + 100)]);
													gameForListeners.placing.alpha = 0.5;
													break;
												case 2:
													gameForListeners.placing = new Rivet(eventPosition.copy());
													gameForListeners.placing.alpha = 0.5;
													break;
												case 3:
													gameForListeners.placing = new Hinge(eventPosition.copy());
													gameForListeners.placing.alpha = 0.5;
													break;
											}
										}
									}
								}
							} else {
								for (var i=gameForListeners.bodies.length-1; i>=0; i--) {
									let hitLink = false;
									for (var j=gameForListeners.bodies[i].rivets.length-1; j>=0; j--) {
										if (gameForListeners.bodies[i].rivets[j].isInside(eventPosition)) {
											gameForListeners.mousePosition.set(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y);
											gameForListeners.dragging = gameForListeners.bodies[i].rivets[j];
											hitLink = true;
											break;
										}
									}

									if (hitLink) {
										break;
									}

									for (var j=gameForListeners.bodies[i].hinges.length-1; j>=0; j--) {
										if (gameForListeners.bodies[i].hinges[j].isInside(eventPosition)) {
											gameForListeners.mousePosition.set(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y);
											gameForListeners.dragging = gameForListeners.bodies[i].hinges[j];
											hitLink = true;
											break;
										}
									}

									if (hitLink) {
										break;
									}

									if (gameForListeners.bodies[i].isInside(eventPosition)) {
										gameForListeners.mousePosition.set(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y);
										gameForListeners.dragging = gameForListeners.bodies[i];
										break;
									}
								}
							}
						}
					}
					break;
				case 3:
					gameForListeners.contextMenu = null;
					if (gameForListeners.canInteract && gameForListeners.dragging == null && gameForListeners.placing == null) {
						for (var i=gameForListeners.bodies.length-1; i>=0; i--) {
							let hitLink = false;
							for (var j=gameForListeners.bodies[i].rivets.length-1; j>=0; j--) {
								if (gameForListeners.bodies[i].rivets[j].isInside(eventPosition)) {
									let contextMenuPosition = eventPosition.copy();
									let distanceFromTop = contextMenuPosition.y - gameForListeners.contextMenuHeight - 2*gameForListeners.contextMenuCurveSize - 2;
									let distanceFromSide = window.innerWidth - (contextMenuPosition.x + gameForListeners.contextMenuWidth + 2*gameForListeners.contextMenuCurveSize + 2);
									if (distanceFromTop < 0) {
										contextMenuPosition.addVector(new Vector(0, -distanceFromTop));
									}
									if (distanceFromSide < 0) {
										contextMenuPosition.addVector(new Vector(distanceFromSide, 0));
									}

									gameForListeners.contextMenu = [gameForListeners.bodies[i].rivets[j], contextMenuPosition];
									hitLink = true;
									break;
								}
							}

							if (hitLink) {
								break;
							}

							for (var j=gameForListeners.bodies[i].hinges.length-1; j>=0; j--) {
								if (gameForListeners.bodies[i].hinges[j].isInside(eventPosition)) {
									let contextMenuPosition = eventPosition.copy();
									let distanceFromTop = contextMenuPosition.y - gameForListeners.contextMenuHeight - 2*gameForListeners.contextMenuCurveSize - 2;
									let distanceFromSide = window.innerWidth - (contextMenuPosition.x + gameForListeners.contextMenuWidth + 2*gameForListeners.contextMenuCurveSize + 2);
									if (distanceFromTop < 0) {
										contextMenuPosition.addVector(new Vector(0, -distanceFromTop));
									}
									if (distanceFromSide < 0) {
										contextMenuPosition.addVector(new Vector(distanceFromSide, 0));
									}

									gameForListeners.contextMenu = [gameForListeners.bodies[i].hinges[j], contextMenuPosition];
									hitLink = true;
									break;
								}
							}

							if (hitLink) {
								break;
							}

							if (gameForListeners.bodies[i].isInside(eventPosition)) {
								let contextMenuPosition = eventPosition.copy();
								let distanceFromTop = contextMenuPosition.y - gameForListeners.contextMenuHeight - 2*gameForListeners.contextMenuCurveSize - 2;
								let distanceFromSide = window.innerWidth - (contextMenuPosition.x + gameForListeners.contextMenuWidth + 2*gameForListeners.contextMenuCurveSize + 2);
								if (distanceFromTop < 0) {
									contextMenuPosition.addVector(new Vector(0, -distanceFromTop));
								}
								if (distanceFromSide < 0) {
									contextMenuPosition.addVector(new Vector(distanceFromSide, 0));
								}

								gameForListeners.contextMenu = [gameForListeners.bodies[i], contextMenuPosition];
								break;
							}
						}
					}
					break;
			}
		});

		addMouseUpListener(function(which, eventX, eventY) {
			switch(which) {
				case 1:
					let eventPosition = new Vector(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y);

					gameForListeners.dragging = null;
					gameForListeners.draggingSlider = null;

					if (gameForListeners.placing != null && gameForListeners.placingDrag) {
						if (gameForListeners.placing instanceof Body && getDistance(new Vector(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y), gameForListeners.placing.vertices[0]) > 5) {
							gameForListeners.placing.computeArea();
							gameForListeners.placing.computeMass();
							if (gameForListeners.placing.vertices.length > 1) {
								gameForListeners.placing.computeNormals();
							}

							gameForListeners.addBody(gameForListeners.placing);

							gameForListeners.placing.alpha = 1;
							gameForListeners.placing = null;
							gameForListeners.placingDrag = false;
						} else if (gameForListeners.placing instanceof Link) {
							gameForListeners.addLink(gameForListeners.placing);

							gameForListeners.placing.alpha = 1;
							gameForListeners.placing = null;
							gameForListeners.placingDrag = false;
						}
					}
					break;
			}
		});

		addMouseMoveListener(function(eventX, eventY) {
			gameForListeners.mousePosition.set(eventX + gameForListeners.camera.x, eventY + gameForListeners.camera.y);
		});

		addKeyDownListener(function(keyCode) {
			if (keyCode != 'Space') {
				if (!contains(gameForListeners.inputs, keyCode)) {
					gameForListeners.inputs.push(keyCode);
				}
			}
			switch(keyCode) {
				case 'Space':
					gameForListeners.changePause();
					break;
			}
		});

		addKeyUpListener(function(keyCode) {
			if (keyCode != 'Space') {
				remove(gameForListeners.inputs, keyCode);
			}
		});
	}

	addBody(body) {
		this.bodies.push(body);
		body.id = this.currid;
		this.currid++;
	}

	removeBody(body) {
		if (this.contextMenu != null && this.contextMenu[0] == body) {
			this.contextMenu = null;
		}
		if (this.dragging != null && this.dragging == body) {
			this.dragging = null;
		}
		if (this.placing != null && this.placing == body) {
			this.placing = null;
			this.placingDrag = false;
		}
		remove(this.bodies, body);
		body.delete();
	}

	addLink(link) {
		link.findLinks(this.bodies);
	}

	removeLink(link) {
		if (this.contextMenu != null && this.contextMenu[0] == link) {
			this.contextMenu = null;
		}
		if (this.dragging != null && this.dragging == link) {
			this.dragging = null;
		}
		if (this.placing != null && this.placing == link) {
			this.placing = null;
			this.placingDrag = false;
		}
		link.delete();
	}

	changePause() {
		this.pause = !this.pause;
		this.canInteract = this.pause;

		if (!this.canInteract) {
			this.dragging = null;
			this.contextMenu = null;
			this.draggingSlider = null;
		}
	}

	setStartTime() {
		this.startTime = new Date();
	}

	timeSinceStart() {
		return new Date().getTime() - this.startTime;
	}

	tick() {
		this.tickID++;

		if (this.canInteract) {
			if (this.placing != null) {
				if (this.placingDrag) {
					if (this.placing instanceof Body) {
						if (this.placing instanceof Circle) {
							this.placing.radius = Math.round(getDistance(this.mousePosition, this.placing.vertices[0]));
						} else if (this.placing instanceof Rectangle) {
							let xDiff = this.mousePosition.x - this.placing.vertices[0].x;
							let yDiff = this.mousePosition.y - this.placing.vertices[0].y;
							if (xDiff * yDiff > 0) {
								this.placing.vertices[1].set(this.mousePosition.x, this.placing.vertices[0].y);
								this.placing.vertices[2].set(this.mousePosition.x, this.mousePosition.y);
								this.placing.vertices[3].set(this.placing.vertices[0].x, this.mousePosition.y);
							} else {
								this.placing.vertices[1].set(this.placing.vertices[0].x, this.mousePosition.y);
								this.placing.vertices[2].set(this.mousePosition.x, this.mousePosition.y);
								this.placing.vertices[3].set(this.mousePosition.x, this.placing.vertices[0].y);
							}
						}
					}
				} else {
					this.placing.translate(this.mousePosition.differenceVector(this.oldMousePosition));
				}
			} else if (this.dragging != null) {
				let movementVector = this.mousePosition.differenceVector(this.oldMousePosition);
				this.dragging.translate(movementVector);

				/*let startBoxMinX = this.startBox.vertices[0].x;
				let startBoxMaxX = this.startBox.vertices[0].x;
				let startBoxMinY = this.startBox.vertices[0].y;
				let startBoxMaxY = this.startBox.vertices[0].y;

				for (var i=0; i<this.startBox.vertices.length; i++) {
					if (this.startBox.vertices[i].x < startBoxMinX) {
						startBoxMinX = this.startBox.vertices[i].x;
					}
					if (this.startBox.vertices[i].x > startBoxMaxX) {
						startBoxMaxX = this.startBox.vertices[i].x;
					}
					if (this.startBox.vertices[i].y < startBoxMinY) {
						startBoxMinY = this.startBox.vertices[i].y;
					}
					if (this.startBox.vertices[i].y > startBoxMinY) {
						startBoxMaxY = this.startBox.vertices[i].y;
					}
				}

				let illegal = false;
				if (this.dragging instanceof Body) {
					for (var i=0; i<this.dragging.vertices.length; i++) {
						if (this.dragging.vertices[i].x < startBoxMinX || 
							this.dragging.vertices[i].x > startBoxMaxX ||
							this.dragging.vertices[i].y < startBoxMinY ||
							this.dragging.vertices[i].y > startBoxMaxY) {
							illegal = true;
							break;
						}
					}
				} else if (this.dragging instanceof Link) {
					if (this.dragging.position.x < startBoxMinX || 
						this.dragging.position.x > startBoxMaxX ||
						this.dragging.position.y < startBoxMinY ||
						this.dragging.position.y > startBoxMaxY) {
						illegal = true;
					}
				}

				if (illegal) {
					movementVector.negate();
					this.dragging.translate(movementVector);
				}*/

				if (this.dragging instanceof Body) {
					let moved = false;
					for (var i=0; i<this.dragging.rivets.length; i++) {
						if (!this.dragging.isInside(this.dragging.rivets[i].position)) {
							movementVector.negate();
							this.dragging.translate(movementVector);
							moved = true;
							break;
						}
					}

					if (!moved) {
						for (var i=0; i<this.dragging.hinges.length; i++) {
							if (!this.dragging.isInside(this.dragging.hinges[i].position)) {
								movementVector.negate();
								this.dragging.translate(movementVector);
								break;
							}
						}
					}
				} else if (this.dragging instanceof Link) {
					if (!(this.dragging.body1 == null || this.dragging.body1.isInside(this.dragging.position)) || !(this.dragging.body2 == null || this.dragging.body2.isInside(this.dragging.position))) {
						movementVector.negate();
						this.dragging.translate(movementVector);
					}
				}
			} else if (this.contextMenu != null && this.draggingSlider != null) {
				let sliders = [];
				let slidersMin = [];
				let slidersMax = [];
				if (this.contextMenu[0] instanceof Body) {
					sliders = this.bodySliders;
					slidersMin = this.bodySlidersMin;
					slidersMax = this.bodySlidersMax;
				} else if (this.contextMenu[0] instanceof Rivet) {
					sliders = this.rivetSliders;
					slidersMin = this.rivetSlidersMin;
					slidersMax = this.rivetSlidersMax;
				} else if (this.contextMenu[0] instanceof Hinge) {
					sliders = this.hingeSliders;
					slidersMin = this.hingeSlidersMin;
					slidersMax = this.hingeSlidersMax;
				}

				let sliderStartX = this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2;
				let percentage = (this.mousePosition.x - sliderStartX)/this.contextMenuSliderWidth;
				let actualNumber = Math.min(Math.max(percentage, 0), 1) * (slidersMax[this.draggingSlider] - slidersMin[this.draggingSlider]) + slidersMin[this.draggingSlider];
				if (Math.abs(slidersMax[this.draggingSlider]) <= 10) {
					actualNumber = Math.round(actualNumber * 1000)/1000;
				} else {
					actualNumber = Math.round(actualNumber);
				}

				if (sliders[this.draggingSlider] == 'angle') {
					this.contextMenu[0].rotate(actualNumber - this.contextMenu[0]['angle']);
				} else if (sliders[this.draggingSlider] == 'density') {
					this.contextMenu[0][sliders[this.draggingSlider]] = actualNumber;
					this.contextMenu[0].computeMass();
				} else {
					this.contextMenu[0][sliders[this.draggingSlider]] = actualNumber;
				}
			}
		}

		if (!this.pause) {
			this.collisions = [];
			let hinges = [];
			for (var i=0; i<this.bodies.length; i++) {
				for (var j=i+1; j<this.bodies.length; j++) {
					let hasHinge = false;
					for (var k=0; k<this.bodies[i].hinges.length; k++) {
						if (this.bodies[i].hinges[k].getPartner(this.bodies[i]) == this.bodies[j]) {
							hasHinge = true;
							break;
						}
					}

					if (hasHinge) {
						break;
					}

					if ((this.bodies[i].compositeObject == null && this.bodies[i].compositeObject == null) || this.bodies[i].compositeObject != this.bodies[j].compositeObject) {
						let collision = new Collision(this.bodies[i], this.bodies[j]);
						collision.solve();

						if (collision.contacts.length > 0) {
							this.collisions.push(collision);
						}
					}
				}

				this.bodies[i].applyHingeImpulses();
			}

			for (var i=0; i<this.collisions.length; i++) {
				this.collisions[i].applyImpulse();
			}

			for (var i=0; i<this.bodies.length; i++) {
				this.bodies[i].tick(this.dt);
			}

			for (var i=0; i<this.collisions.length; i++) {
				this.collisions[i].positionalCorrection();
			}
		}

		this.oldMousePosition.set(this.mousePosition.x, this.mousePosition.y);
	}

	render() {
		for (var i=0; i<this.inputs.length; i++) {
			switch(this.inputs[i]) {
				case 'KeyW':
					this.camera.y += -5;
					break;
				case 'KeyA':
					this.camera.x += -5;
					break;
				case 'KeyS':
					this.camera.y += 5;
					break;
				case 'KeyD':
					this.camera.x += 5;
					break;
				case 'ArrowRight':
					this.camera.x += 5;
					break;
				case 'ArrowLeft':
					this.camera.x += -5;
					break;
				case 'ArrowUp':
					this.camera.y += -5;
					break;
				case 'ArrowDown':
					this.camera.y += 5;
					break;
			}
		}

		this.context.fillStyle = 'rgba(115, 215, 255, 1)';
		this.context.fillRect(0, 0, this.canvas.width, this.canvas.height);

		for (var i=0; i<this.bodies.length; i++) {
			this.bodies[i].render(this.context, this.camera);

			for (var j=0; j<this.bodies[i].rivets.length; j++) {
				this.bodies[i].rivets[j].render(this.context, this.camera);
			}

			for (var j=0; j<this.bodies[i].hinges.length; j++) {
				this.bodies[i].hinges[j].render(this.context, this.camera);
			}
		}

		this.context.fillStyle = 'rgba(40, 40, 40, 0.7)';
		this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
		this.context.lineWidth = 1;
		this.context.beginPath();
		this.context.moveTo(0, this.buildMenuY - this.buildMenuCurveSize);
		this.context.lineTo(this.buildMenuWidth, this.buildMenuY - this.buildMenuCurveSize);
		this.context.arc(this.buildMenuWidth, this.buildMenuY, this.buildMenuCurveSize, 3*Math.PI/2, 2*Math.PI, false);
		this.context.lineTo(this.buildMenuWidth + this.buildMenuCurveSize, this.buildMenuY + this.buildMenuHeight + this.buildMenuCurveSize);
		this.context.arc(this.buildMenuWidth, this.buildMenuY + this.buildMenuHeight + this.buildMenuCurveSize, this.buildMenuCurveSize, 0, Math.PI/2, false);
		this.context.lineTo(0, this.buildMenuY + this.buildMenuHeight + 2*this.buildMenuCurveSize);
		this.context.fill();
		this.context.stroke();
		this.context.closePath();

		let buildMenuIcons = 0;

		this.context.font = 'bold 16px sans-serif';
		this.context.textAlign = 'center';
		this.context.textBaseline = 'alphabetic';
		this.context.fillStyle = 'rgba(0, 0, 0, 1)';
		this.context.fillText('Circle', (this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize - 2);
		this.context.fillStyle = 'rgba(255, 255, 255, 1)';
		this.context.fillRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize, this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.strokeRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize, this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.fillStyle = 'rgba(255, 0, 0, 1)';
		this.context.beginPath();
		this.context.arc((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize + this.buildMenuIconSize/2, this.buildMenuIconSize/3, 0, 2*Math.PI, false);
		this.context.fill();
		this.context.stroke();
		this.context.closePath();
		buildMenuIcons++;

		this.context.font = 'bold 16px sans-serif';
		this.context.textAlign = 'center';
		this.context.textBaseline = 'alphabetic';
		this.context.fillStyle = 'rgba(0, 0, 0, 1)';
		this.context.fillText('Rect', (this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize - 2 + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize));
		this.context.fillStyle = 'rgba(255, 255, 255, 1)';
		this.context.fillRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.strokeRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.fillStyle = 'rgba(255, 0, 0, 1)';
		this.context.fillRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2 + this.buildMenuIconSize/10, this.buildMenuY + this.buildMenuCurveSize + this.buildMenuIconSize/4 + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), 
								this.buildMenuIconSize*0.8, this.buildMenuIconSize*0.5);
		this.context.strokeRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2 + this.buildMenuIconSize/10, this.buildMenuY + this.buildMenuCurveSize + this.buildMenuIconSize/4 + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), 
								this.buildMenuIconSize*0.8, this.buildMenuIconSize*0.5);
		this.context.closePath();
		buildMenuIcons++;

		this.context.font = 'bold 16px sans-serif';
		this.context.textAlign = 'center';
		this.context.textBaseline = 'alphabetic';
		this.context.fillStyle = 'rgba(0, 0, 0, 1)';
		this.context.fillText('Rivet', (this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize - 2 + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize));
		this.context.fillStyle = 'rgba(255, 255, 255, 1)';
		this.context.fillRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.strokeRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.strokeStyle = 'rgba(0, 0, 255, 1)';
		this.context.lineWidth = 4;
		this.context.beginPath();
		this.context.moveTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2 + this.buildMenuIconSize*0.2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.2);
		this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2 + this.buildMenuIconSize*0.8, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.8);
		this.context.moveTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2 + this.buildMenuIconSize*0.8, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.2);
		this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2 + this.buildMenuIconSize*0.2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.8);
		this.context.stroke();
		this.context.closePath();
		buildMenuIcons++;

		this.context.font = 'bold 16px sans-serif';
		this.context.textAlign = 'center';
		this.context.textBaseline = 'alphabetic';
		this.context.fillStyle = 'rgba(0, 0, 0, 1)';
		this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
		this.context.lineWidth = 1;
		this.context.fillText('Hinge', (this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize - 2 + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize));
		this.context.fillStyle = 'rgba(255, 255, 255, 1)';
		this.context.fillRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.strokeRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
		this.context.strokeStyle = 'rgba(0, 0, 255, 1)';
		this.context.lineWidth = 4;
		this.context.beginPath();
		this.context.arc((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.5,
							this.buildMenuIconSize * 0.4, 0, 2*Math.PI, false);
		this.context.stroke();
		this.context.closePath();
		buildMenuIcons++;

		if (this.contextMenu != null) {
			let top = this.contextMenu[1].y - 2*this.contextMenuCurveSize - this.contextMenuHeight;

			this.context.fillStyle = 'rgba(40, 40, 40, 0.7)';
			this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
			this.context.lineWidth = 1;
			this.context.beginPath();
			this.context.moveTo(this.contextMenu[1].x + this.contextMenuCurveSize, this.contextMenu[1].y - 2*this.contextMenuCurveSize - this.contextMenuHeight);
			this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + this.contextMenuWidth, this.contextMenu[1].y - 2*this.contextMenuCurveSize - this.contextMenuHeight);
			this.context.arc(this.contextMenu[1].x + this.contextMenuCurveSize + this.contextMenuWidth, this.contextMenu[1].y - this.contextMenuCurveSize - this.contextMenuHeight, this.contextMenuCurveSize, 3*Math.PI/2, 2*Math.PI, false);
			this.context.lineTo(this.contextMenu[1].x + 2*this.contextMenuCurveSize + this.contextMenuWidth, this.contextMenu[1].y - this.contextMenuCurveSize);
			this.context.arc(this.contextMenu[1].x + this.contextMenuCurveSize + this.contextMenuWidth, this.contextMenu[1].y - this.contextMenuCurveSize, this.contextMenuCurveSize, 0, Math.PI/2, false);
			this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize, this.contextMenu[1].y);
			this.context.arc(this.contextMenu[1].x + this.contextMenuCurveSize, this.contextMenu[1].y - this.contextMenuCurveSize, this.contextMenuCurveSize, Math.PI/2, Math.PI, false);
			this.context.lineTo(this.contextMenu[1].x, this.contextMenu[1].y - this.contextMenuCurveSize - this.contextMenuHeight);
			this.context.arc(this.contextMenu[1].x + this.contextMenuCurveSize, this.contextMenu[1].y - this.contextMenuCurveSize - this.contextMenuHeight, this.contextMenuCurveSize, Math.PI, 3*Math.PI/2, false);
			this.context.fill();
			this.context.stroke();
			this.context.closePath();

			this.context.font = 'bold 28px sans-serif';
			this.context.textAlign = 'left';
			this.context.textBaseline = 'top';
			this.context.fillStyle = 'rgba(190, 0, 0, 1)';
			this.context.fillRect(this.contextMenu[1].x + this.contextMenuCurveSize, top + this.contextMenuCurveSize + 22, this.contextMenuWidth, 3);

			this.context.fillStyle = 'rgba(255, 50, 50, 1)';
			this.context.lineWidth = 1;
			this.context.fillText(this.contextMenu[0].getName(), this.contextMenu[1].x + this.contextMenuCurveSize, top + this.contextMenuCurveSize - 2);

			this.context.strokeStyle = 'rgba(255, 0, 0, 1)';
			this.context.lineWidth = 4;
			this.context.beginPath();
			this.context.moveTo(this.contextMenu[1].x + 2*this.contextMenuCurveSize + this.contextMenuWidth - 12, top + 8);
			this.context.lineTo(this.contextMenu[1].x + 2*this.contextMenuCurveSize + this.contextMenuWidth - 12 - this.contextMenuCurveSize, top + 8 + this.contextMenuCurveSize);
			this.context.moveTo(this.contextMenu[1].x + 2*this.contextMenuCurveSize + this.contextMenuWidth - 12 - this.contextMenuCurveSize, top + 8);
			this.context.lineTo(this.contextMenu[1].x + 2*this.contextMenuCurveSize + this.contextMenuWidth - 12, top + 8 + this.contextMenuCurveSize);
			this.context.stroke();
			this.context.closePath();
			this.context.fillStyle = 'rgba(0, 0, 0, 1)';
			this.context.font = '12px sans-serif';
			this.context.textAlign = 'center';
			this.context.textBaseline = 'top';
			this.context.fillText('Delete', this.contextMenu[1].x + 2*this.contextMenuCurveSize + this.contextMenuWidth - 12 - this.contextMenuCurveSize/2, top + 8 + this.contextMenuCurveSize + 2);

			let sliders = [];
			let slidersMin = [];
			let slidersMax = [];
			if (this.contextMenu[0] instanceof Body) {
				sliders = this.bodySliders;
				slidersMin = this.bodySlidersMin;
				slidersMax = this.bodySlidersMax;
			} else if (this.contextMenu[0] instanceof Rivet) {
				sliders = this.rivetSliders;
				slidersMin = this.rivetSlidersMin;
				slidersMax = this.rivetSlidersMax;

				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = 'bold 16px sans-serif';
				this.context.textAlign = 'center';
				this.context.textBaseline = 'top';
				this.context.fillText('Connects: ' + (this.contextMenu[0].body1 == null ? 'Background' : this.contextMenu[0].body1.getName()) + ' to ' + (this.contextMenu[0].body2 == null ? 'Background' : this.contextMenu[0].body2.getName()),
										this.contextMenu[1].x + this.contextMenuCurveSize + this.contextMenuWidth/2, this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length));
			} else if (this.contextMenu[0] instanceof Hinge) {
				sliders = this.hingeSliders;
				slidersMin = this.hingeSlidersMin;
				slidersMax = this.hingeSlidersMax;

				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = 'bold 16px sans-serif';
				this.context.textAlign = 'center';
				this.context.textBaseline = 'top';
				this.context.fillText('Connects: ' + (this.contextMenu[0].body1 == null ? 'Background' : this.contextMenu[0].body1.getName()) + ' to ' + (this.contextMenu[0].body2 == null ? 'Background' : this.contextMenu[0].body2.getName()),
										this.contextMenu[1].x + this.contextMenuCurveSize + this.contextMenuWidth/2, this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length));
			}

			let sliderHeight = 0;
			for (var j=0; j<sliders.length; j++) {
				let variable = this.contextMenu[0][sliders[j]];
				sliderHeight = this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * j;
				this.context.font = '12px sans-serif';
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
				this.context.textAlign = 'center';
				this.context.textBaseline = 'alphabetic';
				this.context.fillText(sliders[j][0].toUpperCase() + sliders[j].slice(1), this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2, this.contextMenu[1].y - sliderHeight - 10);
				this.context.textBaseline = 'top';
				this.context.fillText(slidersMin[j], this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2, this.contextMenu[1].y - sliderHeight + 8);
				this.context.fillText(slidersMax[j], this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + this.contextMenuSliderWidth, this.contextMenu[1].y - sliderHeight + 8);
				this.context.fillStyle = 'rgba(220, 220, 220, 1)';
				this.context.fillText(variable, this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + this.contextMenuSliderWidth/2, this.contextMenu[1].y - sliderHeight + 8);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillRect(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2, this.contextMenu[1].y - sliderHeight, this.contextMenuSliderWidth, this.contextMenuSliderThickness);
				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.lineWidth = 1;
				this.context.beginPath();
				this.context.moveTo(Math.round(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + ((variable - slidersMin[j])/(slidersMax[j] - slidersMin[j])) * this.contextMenuSliderWidth - this.contextMenuSliderGrabberWidth/2), 
									Math.round(this.contextMenu[1].y - sliderHeight + 2 - this.contextMenuSliderGrabberHeight/2));
				this.context.lineTo(Math.round(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + ((variable - slidersMin[j])/(slidersMax[j] - slidersMin[j])) * this.contextMenuSliderWidth + this.contextMenuSliderGrabberWidth/2), 
									Math.round(this.contextMenu[1].y - sliderHeight + 2 - this.contextMenuSliderGrabberHeight/2));
				this.context.lineTo(Math.round(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + ((variable - slidersMin[j])/(slidersMax[j] - slidersMin[j])) * this.contextMenuSliderWidth + this.contextMenuSliderGrabberWidth/2), 
									Math.round(this.contextMenu[1].y - sliderHeight + 2));
				this.context.lineTo(Math.round(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + ((variable - slidersMin[j])/(slidersMax[j] - slidersMin[j])) * this.contextMenuSliderWidth), 
									Math.round(this.contextMenu[1].y - sliderHeight + 2 + this.contextMenuSliderGrabberHeight/2));
				this.context.lineTo(Math.round(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + ((variable - slidersMin[j])/(slidersMax[j] - slidersMin[j])) * this.contextMenuSliderWidth - this.contextMenuSliderGrabberWidth/2), 
									Math.round(this.contextMenu[1].y - sliderHeight + 2));
				this.context.lineTo(Math.round(this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 + ((variable - slidersMin[j])/(slidersMax[j] - slidersMin[j])) * this.contextMenuSliderWidth - this.contextMenuSliderGrabberWidth/2), 
									Math.round(this.contextMenu[1].y - sliderHeight + 2 - this.contextMenuSliderGrabberHeight/2));
				this.context.fill();
				this.context.stroke();
				this.context.closePath();
			}
		}

		if (this.placing != null) {
			this.placing.render(this.context, this.camera);
		}

		if (this.debug) {
			for (var i=0; i<this.collisions.length; i++) {
				for (var j=0; j<this.collisions[i].contacts.length; j++) {
					this.context.strokeStyle = 'rgba(0, 255, 0, 1)';
					this.context.lineWidth = 5;
					this.context.beginPath();
					this.context.arc(this.collisions[i].contacts[j].x, this.collisions[i].contacts[j].y, 10, 0, 2*Math.PI, false);
					this.context.stroke();
					this.context.closePath();
				}
			}
		}
	}
}

let totalTimeUnloaded = 0;
let mostRecentUnload = null;
document.addEventListener('visibilitychange', function () {
	if (document.visibilityState === 'visible' && mostRecentUnload) {
		totalTimeUnloaded += (new Date().getTime() - mostRecentUnload.getTime());
		mostRecentUnload = null;
	} else if (document.visibilityState === 'hidden' && !mostRecentUnload) {
		mostRecentUnload = new Date();
	}
});

function gameLoop(game) {
	while (game.timeSinceStart() - totalTimeUnloaded >= (1000/game.ticksPerSecond) * game.tickID) {
		game.tick();
	}
	game.render();

	window.requestAnimationFrame(function() {gameLoop(game);});
}

function start() {
	let game = new Game();
	game.setStartTime();
	window.requestAnimationFrame(function() {gameLoop(game);});
	window.onload = function() {
		document.body.appendChild(game.canvas);
	}
}

start();