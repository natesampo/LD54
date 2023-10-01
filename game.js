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

		this.body1Index = -1;
		this.body2Index = -1;

		this.partOfLevel = false;
	}

	place() {
		this.alpha = 1;
	}

	findLinks(bodies) {
		if (this.body1 == null && this.body2 == null) {
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

	delete() {
	}

	exportString(index) {
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

		newContext.strokeStyle = 'rgba(0, 0, 0, 1)';
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

	exportString(index) {
		let shapeExportString = ',l,' + this.position.x + ',' + this.position.y + ',' + this.body1Index + ',' + this.body2Index + ',' + this.size + ',' + this.angle.toFixed(3) + ',' + this.red + ',' + this.green + ',' + this.blue + ',' + this.alpha;

		if (this.partOfLevel) {
			shapeExportString += ',1';
		} else {
			shapeExportString += ',0';
		}

		return shapeExportString;
	}
}

class Spring extends Link {
	constructor(position) {
		super(position);
		this.blue = 255;
		this.strength = 10;
		this.targetLength = 50;
		this.damping = 1;
		this.width = 10;

		this.position2 = new Vector(this.position.x, this.position.y - 50);

		this.body1Distance = null;
		this.body2Distance = null;
		this.body1Theta = null;
		this.body2Theta = null;
		this.body1InitialAngle = null;
		this.body2InitialAngle = null;

		this.lastDisplacement = null;
	}

	place() {
		super.place();

		this.targetLength = Math.ceil(this.position2.differenceVector(this.position).magnitude());
		if (this.targetLength == 0) {
			this.targetLength = 1;
		}
	}

	translate(vector) {
		super.translate(vector);

		this.position2.translate(vector);
	}

	render(context, camera) {
		let canvas = context.canvas;

		context.fillStyle = 'rgba(0, 0, 0, 1)';
		context.beginPath();
		context.arc(this.position.x - camera.x, this.position.y - camera.y, 3, 0, 2*Math.PI);
		context.fill();
		context.closePath();

		context.beginPath();
		context.arc(this.position2.x - camera.x, this.position2.y - camera.y, 3, 0, 2*Math.PI);
		context.fill();
		context.closePath();

		let numLoops = Math.ceil(this.targetLength/15);
		if (numLoops > 0) {
			context.lineWidth = 4;
			let distanceBetweenLoops = this.position2.differenceVector(this.position).magnitude() / numLoops;
			let angle = getAngle(this.position, this.position2);
			let cos = Math.cos(angle);
			let sin = Math.sin(angle);
			let cos2 = Math.cos(angle + Math.PI/2);
			let sin2 = Math.sin(angle + Math.PI/2);
			context.strokeStyle = 'rgba(0, 0, 0, ' + this.alpha + ')';
			context.beginPath();
			context.moveTo(this.position.x - camera.x, this.position.y - camera.y);
			for (var i=0; i<numLoops; i++) {
				context.lineTo(this.position.x - camera.x + distanceBetweenLoops * (i+0.5) * cos - this.width * cos2, this.position.y - camera.y + distanceBetweenLoops * (i+0.5) * sin - this.width*sin2);
				context.lineTo(this.position.x - camera.x + distanceBetweenLoops * (i+0.5) * cos + this.width * cos2, this.position.y - camera.y + distanceBetweenLoops * (i+0.5) * sin + this.width*sin2);
			}
			context.lineTo(this.position2.x - camera.x, this.position2.y - camera.y);
			context.stroke();
			context.closePath();

			context.lineWidth = 3;
			context.strokeStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + this.alpha + ')';
			context.beginPath();
			context.moveTo(this.position.x - camera.x, this.position.y - camera.y);
			for (var i=0; i<numLoops; i++) {
				context.lineTo(this.position.x - camera.x + distanceBetweenLoops * (i+0.5) * cos - this.width * cos2, this.position.y - camera.y + distanceBetweenLoops * (i+0.5) * sin - this.width*sin2);
				context.lineTo(this.position.x - camera.x + distanceBetweenLoops * (i+0.5) * cos + this.width * cos2, this.position.y - camera.y + distanceBetweenLoops * (i+0.5) * sin + this.width*sin2);
			}
			context.lineTo(this.position2.x - camera.x, this.position2.y - camera.y);
			context.stroke();
			context.closePath();
		}
	}

	initPoints() {
		if (this.body1 != null) {
			let body1COM = this.body1.findCenterOfMass();
			this.body1Distance = getDistance(body1COM, this.position);
			this.body1Theta = getAngle(body1COM, this.position);
			this.body1InitialAngle = this.body1.angle;
		}

		if (this.body2 != null) {
			let body2COM = this.body2.findCenterOfMass();
			this.body2Distance = getDistance(body2COM, this.position2);
			this.body2Theta = getAngle(body2COM, this.position2);
			this.body2InitialAngle = this.body2.angle;
		}
	}

	applyImpulse() {
		let firstAnchorPoint = this.position.copy();
		if (this.body1 != null) {
			let body1COM = this.body1.findCenterOfMass();
			firstAnchorPoint = new Vector(body1COM.x + this.body1Distance * Math.cos(this.body1.angle - this.body1InitialAngle + this.body1Theta), body1COM.y + this.body1Distance * Math.sin(this.body1.angle - this.body1InitialAngle + this.body1Theta));
		}

		let secondAnchorPoint = this.position2.copy();
		if (this.body2 != null) {
			let body2COM = this.body2.findCenterOfMass();
			secondAnchorPoint = new Vector(body2COM.x + this.body2Distance * Math.cos(this.body2.angle - this.body2InitialAngle + this.body2Theta), body2COM.y + this.body2Distance * Math.sin(this.body2.angle - this.body2InitialAngle + this.body2Theta));
		}

		let anchorPointDiffX = secondAnchorPoint.x - firstAnchorPoint.x;
		let anchorPointDiffY = secondAnchorPoint.y - firstAnchorPoint.y;

		this.position.set(firstAnchorPoint.x, firstAnchorPoint.y);
		this.position2.set(secondAnchorPoint.x, secondAnchorPoint.y);

		let directionVector = firstAnchorPoint.differenceVector(secondAnchorPoint);
		let separation = directionVector.magnitude();
		let displacement = separation - this.targetLength;

		let relativeVelocity = 0;
		if (this.lastDisplacement != null) {
			relativeVelocity = displacement - this.lastDisplacement;
		}

		directionVector.normalize();
		directionVector.multiplyScalar(-this.strength * displacement - relativeVelocity * 50*this.damping);

		if (this.body1 != null) {
			firstAnchorPoint.subtractVector(this.body1.findCenterOfMass());
			this.body1.applyImpulse(directionVector, firstAnchorPoint);
		}
		
		if (this.body2 != null) {
			secondAnchorPoint.subtractVector(this.body2.findCenterOfMass());
			directionVector.negate();
			this.body2.applyImpulse(directionVector, secondAnchorPoint);
		}

		this.lastDisplacement = displacement;
	}

	isInside(vector) {
		if (pointToLineDistance(vector, this.position, this.position2) < this.width + 3) {
			let distance = getDistance(this.position, this.position2);
			if (getDistance(vector, this.position) > distance + 5) {
				return false;
			}
			if (getDistance(vector, this.position2) > distance + 5) {
				return false;
			}

			return true;
		}
	}

	findLinks(bodies) {
		if (this.body1 == null && this.body2 == null) {
			for (var i=bodies.length-1; i>=0; i--) {
				if (bodies[i].isInside(this.position)) {
					if (this.body1 == null) {
						this.body1 = bodies[i];

						if (this.body2 != null) {
							break;
						}
					}
				} else if (bodies[i].isInside(this.position2)) {
					if (this.body2 == null) {
						this.body2 = bodies[i];

						if (this.body1 != null) {
							break;
						}
					}
				}
			}
		}

		if (this.body1 != null) {
			this.body1.springs.push(this);
		}
		if (this.body2 != null) {
			this.body2.springs.push(this);
		}
	}

	delete() {
		if (this.body1 != null) {
			remove(this.body1.springs, this);
		}
		if (this.body2 != null) {
			remove(this.body2.springs, this);
		}

		this.body1 = null;
		this.body2 = null;
	}

	getName() {
		return 'Spring';
	}

	exportString(index) {
		let shapeExportString = ',s,' + this.position.x + ',' + this.position.y + ',' + this.position2.x + ',' + this.position2.y + ','  + this.body1Index + ',' + this.body2Index + ',' + this.strength + ',' + this.damping.toFixed(3) +
			',' + this.width + ',' + this.targetLength + ',' + this.red + ',' + this.green + ',' + this.blue;

		if (this.partOfLevel) {
			shapeExportString += ',1';
		} else {
			shapeExportString += ',0';
		}

		return shapeExportString;
	}
}

class Hinge extends Link {
	constructor(position) {
		super(position);
		this.blue = 255;
		this.strength = 100;
		this.motorSpeed = 0;
		this.damping = 2;

		this.body1Distance = null;
		this.body2Distance = null;
		this.body1Theta = null;
		this.body2Theta = null;
		this.body1InitialAngle = null;
		this.body2InitialAngle = null;

		this.lastDisplacement = null;
	}

	render(context, camera) {
		let canvas = context.canvas;

		context.strokeStyle = 'rgba(0, 0, 0, ' + this.alpha + ')';
		context.lineWidth = 5;
		context.beginPath();
		context.arc(this.position.x - camera.x, this.position.y - camera.y, this.size/2, 0, 2*Math.PI, false);
		context.stroke();
		context.closePath();

		context.strokeStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + this.alpha + ')';
		context.lineWidth = 3;
		context.beginPath();
		context.arc(this.position.x - camera.x, this.position.y - camera.y, this.size/2, 0, 2*Math.PI, false);
		context.stroke();
		context.closePath();
	}

	initPoints() {
		let body1COM = this.body1.findCenterOfMass();
		this.body1Distance = getDistance(body1COM, this.position);
		this.body1Theta = getAngle(body1COM, this.position);
		this.body1InitialAngle = this.body1.angle;

		if (this.body2 != null) {
			let body2COM = this.body2.findCenterOfMass();
			this.body2Distance = getDistance(body2COM, this.position);
			this.body2Theta = getAngle(body2COM, this.position);
			this.body2InitialAngle = this.body2.angle;
		}
	}

	applyImpulse() {
		let body1COM = this.body1.findCenterOfMass();
		let firstAnchorPoint = new Vector(body1COM.x + this.body1Distance * Math.cos(this.body1.angle - this.body1InitialAngle + this.body1Theta), body1COM.y + this.body1Distance * Math.sin(this.body1.angle - this.body1InitialAngle + this.body1Theta));

		let secondAnchorPoint = this.position.copy();
		if (this.body2 != null) {
			let body2COM = this.body2.findCenterOfMass();
			secondAnchorPoint = new Vector(body2COM.x + this.body2Distance * Math.cos(this.body2.angle - this.body2InitialAngle + this.body2Theta), body2COM.y + this.body2Distance * Math.sin(this.body2.angle - this.body2InitialAngle + this.body2Theta));
		}

		let anchorPointDiffX = secondAnchorPoint.x - firstAnchorPoint.x;
		let anchorPointDiffY = secondAnchorPoint.y - firstAnchorPoint.y;

		if (this.body2 != null) {
			if (this.body1.getMass() > 0 && this.body2.getMass() > 0) {
				let massRatio = this.body1.getMass()/(this.body1.getMass() + this.body2.getMass());
				this.position.set(firstAnchorPoint.x + anchorPointDiffX * massRatio, firstAnchorPoint.y + anchorPointDiffY * massRatio);
			} else if (this.body1.getMass() > 0) {
				this.position.set(secondAnchorPoint.x, secondAnchorPoint.y);
			} else if (this.body2.getMass() > 0) {
				this.position.set(firstAnchorPoint.x, firstAnchorPoint.y);
			} else {
				this.position.set(firstAnchorPoint.x + anchorPointDiffX/2, firstAnchorPoint.y + anchorPointDiffY/2);
			}
		}

		let directionVector = firstAnchorPoint.differenceVector(secondAnchorPoint);
		let displacement = directionVector.magnitude();

		let relativeVelocity = 0;
		if (this.lastDisplacement != null) {
			relativeVelocity = displacement - this.lastDisplacement;
		}

		directionVector.normalize();
		directionVector.multiplyScalar(-this.strength * displacement - relativeVelocity * 50*this.damping);

		firstAnchorPoint.subtractVector(body1COM);
		this.body1.applyImpulse(directionVector, firstAnchorPoint);

		if (this.body2 != null) {
			secondAnchorPoint.subtractVector(this.body2.findCenterOfMass());
			directionVector.negate();
			this.body2.applyImpulse(directionVector, secondAnchorPoint);
		}
		
		if (this.motorSpeed != 0) {
			this.spinMeRightRound();
		}

		this.lastDisplacement = displacement;
	}

	spinMeRightRound() {
		let speedIncrease = 500*this.motorSpeed;
		let maxSpeedDiff = this.motorSpeed/75;

		if (this.body1 != null && this.body2 != null) {
			if ((maxSpeedDiff > 0 && this.body1.angularVelocity - this.body2.angularVelocity < maxSpeedDiff) || (maxSpeedDiff < 0 && this.body1.angularVelocity - this.body2.angularVelocity > maxSpeedDiff)) {
				if (this.body1.getMass() == 0) {
					this.body2.torque += speedIncrease;
				} else if (this.body2.getMass() == 0) {
					this.body1.torque += speedIncrease;
				} else {
					this.body1.torque += speedIncrease * (1-(this.body1.getMass()/(this.body2.getMass() + this.body1.getMass())));
					this.body2.torque -= speedIncrease * (1-(this.body2.getMass()/(this.body2.getMass() + this.body1.getMass())));
				}
			}
		} else {
			if (this.body1 != null && this.body1.getMass() > 0 && ((maxSpeedDiff > 0 && this.body1.angularVelocity < maxSpeedDiff/2) || (maxSpeedDiff < 0 && this.body1.angularVelocity > maxSpeedDiff/2))) {
				this.body1.torque += speedIncrease;
			} else if (this.body2 != null && this.body2.getMass() > 0 && ((maxSpeedDiff > 0 && this.body2.angularVelocity < maxSpeedDiff/2) || (maxSpeedDiff < 0 && this.body2.angularVelocity > maxSpeedDiff/2))) {
				this.body2.torque -= speedIncrease;
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

	exportString(index) {
		let shapeExportString = ',h,' + this.position.x + ',' + this.position.y + ',' + this.body1Index + ',' + this.body2Index + ',' + this.motorSpeed + ',' + this.size + ',' + this.red + ',' + this.green + ',' + this.blue;

		if (this.partOfLevel) {
			shapeExportString += ',1';
		} else {
			shapeExportString += ',0';
		}

		return shapeExportString;
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
		this.friction = 0.25;
		this.red = 255;
		this.green = 0;
		this.blue = 0;
		this.alpha = 1;
		this.rivets = [];
		this.compositeObject = null;
		this.hinges = [];
		this.springs = [];
		this.static = false;
		this.partOfLevel = false;
		this.goal = false;
		this.collisionLayers = {0: true, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};

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

	place() {
		this.alpha = 1;
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
			if (this.hinges[i].body1 == this || this.hinges[i].body1 == null) {
				if ((this.hinges[i].body1 == this && this.hinges[i].body1Distance == null) || (this.hinges[i].body2 == this && this.hinges[i].body2Distance == null)) {
					this.hinges[i].initPoints();
				}

				this.hinges[i].applyImpulse();
			}
		}
	}

	applySpringImpulses() {
		for (var i=0; i<this.springs.length; i++) {
			if (this.springs[i].body1 == this || this.springs[i].body1 == null) {
				if ((this.springs[i].body1 == this && this.springs[i].body1Distance == null) || (this.springs[i].body2 == this && this.springs[i].body2Distance == null)) {
					this.springs[i].initPoints();
				}

				this.springs[i].applyImpulse();
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

	applyTorque(torque) {
		if (this.compositeObject == null) {
			this.angularVelocity += torque;
		} else {
			this.compositeObject.applyTorque(torque);
		}
	}

	setStatic() {
		this.mass = 0;
		this.invMass = 0;
		this.momentOfInertia = 0;
		this.invMomentOfInertia = 0;
		this.static = true;
	}

	isInside(vector) {
		return false;
	}

	delete() {
		for (var i=this.rivets.length-1; i>=0; i--) {
			this.rivets[i].delete();
		}

		for (var i=this.hinges.length-1; i>=0; i--) {
			this.hinges[i].delete();
		}

		for (var i=this.springs.length-1; i>=0; i--) {
			this.springs[i].delete();
		}

		this.rivets = [];
		this.hinges = [];
		this.springs = [];
	}

	tick(dt) {
		if (this.getMass() > 0) {
			let dtForce = this.force.copy();
			dtForce.multiplyScalar(this.getInvMass());
			dtForce.addVector(gravity);
			dtForce.multiplyScalar(dt);
			this.velocity.addVector(dtForce);
			this.translate(this.velocity);

			this.applyTorque(this.torque * dt * this.getInvMass());
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
			this.applyTorque(this.torque * dt * this.getInvMass());
			if (this.angularVelocity != 0) {
				this.rotate(this.angularVelocity);
			}
		}

		this.force.set(0, 0);
		this.torque = 0;
	}

	render(context, camera) {
		let canvas = context.canvas;

		if (this.goal) {
			context.strokeStyle = 'rgba(225, 175, 0, 1)';
			context.fillStyle = 'rgba(255, 215, 0, 1)';
		} else {
			context.strokeStyle = 'rgba(0, 0, 0, 1)';
			context.fillStyle = 'rgba(' + this.red + ', ' + this.green + ', ' + this.blue + ', ' + this.alpha + ')';
		}

		context.lineWidth = 3;
		context.beginPath();
		context.moveTo(this.vertices[0].x - camera.x, this.vertices[0].y - camera.y);
		for (var i=1; i<this.vertices.length; i++) {
			context.lineTo(this.vertices[i].x - camera.x, this.vertices[i].y - camera.y);
		}
		context.lineTo(this.vertices[0].x - camera.x, this.vertices[0].y - camera.y);
		context.lineTo(this.vertices[1].x - camera.x, this.vertices[1].y - camera.y);
		context.fill();
		context.stroke();
		context.closePath();
	}

	getName() {
		return 'Body';
	}

	getShapeExportString(index) {
	}

	exportString(index) {
	}
}

class Circle extends Body {
	constructor(center, radius) {
		super([center]);
		this.radius = radius;

		this.momentOfInertia = 4000000;
		this.invMomentOfInertia = 1/this.momentOfInertia;

		this.computeArea();
		this.computeMass();
	}

	computeArea() {
		if (this.radius) {
			this.area = Math.PI * this.radius * this.radius;
		}
	}

	computeMass() {
		super.computeMass();

		if (this.getMass() > 0) {
			this.momentOfInertia = 4000000 * (this.getMass() / 8000) * (this.getMass() / 8000);
			this.invMomentOfInertia = 1/this.momentOfInertia;
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

	getShapeExportString(index) {
		return ',c,' + this.vertices[0].x.toFixed(3) + ',' + this.vertices[0].y.toFixed(3) + ',' + this.radius.toFixed(3);
	}

	exportString(index) {
		let shapeExportString = this.getShapeExportString(index) + ',' + this.angle.toFixed(3) + ',' + this.density.toFixed(3) + ',' + this.restitution.toFixed(3) + 
			',' + this.friction.toFixed(3) + ',' + this.red + ',' + this.green + ',' + this.blue + ',' + this.alpha;

		let collisionLayerString = ',';
		for (var i in this.collisionLayers) {
			if (this.collisionLayers[i]) {
				collisionLayerString += i;
			}
		}

		let linkExportString = '';
		for (var i=0; i<this.rivets.length; i++) {
			if (this.rivets[i].body1 == this) {
				this.rivets[i].body1Index = index;
				if (this.rivets[i].body2 == null) {
					linkExportString += this.rivets[i].exportString(index);
				} else if (this.rivets[i].body2Index != -1) {
					linkExportString += this.rivets[i].exportString(index);
				}
			} else if (this.rivets[i].body2 == this) {
				this.rivets[i].body2Index = index;
				if (this.rivets[i].body1 == null) {
					linkExportString += this.rivets[i].exportString(index);
				} else if (this.rivets[i].body1Index != -1) {
					linkExportString += this.rivets[i].exportString(index);
				}
			}
		}

		for (var i=0; i<this.hinges.length; i++) {
			if (this.hinges[i].body1 == this) {
				this.hinges[i].body1Index = index;
				if (this.hinges[i].body2 == null) {
					linkExportString += this.hinges[i].exportString(index);
				} else if (this.hinges[i].body2Index != -1) {
					linkExportString += this.hinges[i].exportString(index);
				}
			} else if (this.hinges[i].body2 == this) {
				this.hinges[i].body2Index = index;
				if (this.hinges[i].body1 == null) {
					linkExportString += this.hinges[i].exportString(index);
				} else if (this.hinges[i].body1Index != -1) {
					linkExportString += this.hinges[i].exportString(index);
				}
			}
		}

		for (var i=0; i<this.springs.length; i++) {
			if (this.springs[i].body1 == this) {
				this.springs[i].body1Index = index;
				if (this.springs[i].body2 == null) {
					linkExportString += this.springs[i].exportString(index);
				} else if (this.springs[i].body2Index != -1) {
					linkExportString += this.springs[i].exportString(index);
				}
			} else if (this.springs[i].body2 == this) {
				this.springs[i].body2Index = index;
				if (this.springs[i].body1 == null) {
					linkExportString += this.springs[i].exportString(index);
				} else if (this.springs[i].body1Index != -1) {
					linkExportString += this.springs[i].exportString(index);
				}
			}
		}
		
		let partOfLevelString = ',0';
		if (this.partOfLevel) {
			partOfLevelString = ',1';
		}
		
		return shapeExportString + collisionLayerString + partOfLevelString + linkExportString;
	}
}

class Rectangle extends Body {
	constructor(vertices) {
		super(vertices);
		this.momentOfInertia = 100000000;
		this.invMomentOfInertia = 1/this.momentOfInertia;

		this.computeArea();
		this.computeMass();
	}

	computeArea() {
		this.area = getDistance(this.vertices[0], this.vertices[1]) * getDistance(this.vertices[1], this.vertices[2]);
	}

	computeMass() {
		super.computeMass();

		if (this.getMass() > 0) {
			this.momentOfInertia = 100000000 * (this.getMass() / 20000) * (this.getMass() / 20000);
			this.invMomentOfInertia = 1/this.momentOfInertia;
		}
	}

	isInside(vector) {
		let tempArea = areaOfTriangle(this.vertices[0], vector, this.vertices[3]) + areaOfTriangle(this.vertices[3], vector, this.vertices[2]) + 
						areaOfTriangle(this.vertices[2], vector, this.vertices[1]) + areaOfTriangle(vector, this.vertices[1], this.vertices[0]);

		return Math.floor(tempArea) <= Math.floor(this.area) || Math.abs(Math.floor(tempArea) - Math.floor(this.area)) <= 1;
	}

	getName() {
		return 'Rectangle';
	}

	getShapeExportString(index) {
		let str = this.vertices[0].x.toFixed(3) + ',' + this.vertices[0].y.toFixed(3) + ',' + 
			this.vertices[1].x.toFixed(3) + ',' + this.vertices[1].y.toFixed(3) + ',' + 
			this.vertices[2].x.toFixed(3) + ',' + this.vertices[2].y.toFixed(3) + ',' + 
			this.vertices[3].x.toFixed(3) + ',' + this.vertices[3].y.toFixed(3);

		if (index == -1) {
			return ',st,' + str;
		} else if (index == -2) {
			return ',g,' + str;
		} else if (index == -3) {
			return ',b,' + str;
		}

		return ',r,' + str;
	}

	exportString(index) {
		let shapeExportString = this.getShapeExportString(index) + ',' +
			this.angle.toFixed(3) + ',' + this.density.toFixed(3) + ',' + this.restitution.toFixed(3) + 
			',' + this.friction.toFixed(3) + ',' + this.red + ',' + this.green + ',' + this.blue + ',' + this.alpha;

		let collisionLayerString = ',';
		for (var i in this.collisionLayers) {
			if (this.collisionLayers[i]) {
				collisionLayerString += i;
			}
		}

		let linkExportString = '';
		for (var i=0; i<this.rivets.length; i++) {
			if (this.rivets[i].body1 == this) {
				this.rivets[i].body1Index = index;
				if (this.rivets[i].body2 == null) {
					linkExportString += this.rivets[i].exportString(index);
				} else if (this.rivets[i].body2Index != -1) {
					linkExportString += this.rivets[i].exportString(index);
				}
			} else if (this.rivets[i].body2 == this) {
				this.rivets[i].body2Index = index;
				if (this.rivets[i].body1 == null) {
					linkExportString += this.rivets[i].exportString(index);
				} else if (this.rivets[i].body1Index != -1) {
					linkExportString += this.rivets[i].exportString(index);
				}
			}
		}

		for (var i=0; i<this.hinges.length; i++) {
			if (this.hinges[i].body1 == this) {
				this.hinges[i].body1Index = index;
				if (this.hinges[i].body2 == null) {
					linkExportString += this.hinges[i].exportString(index);
				} else if (this.hinges[i].body2Index != -1) {
					linkExportString += this.hinges[i].exportString(index);
				}
			} else if (this.hinges[i].body2 == this) {
				this.hinges[i].body2Index = index;
				if (this.hinges[i].body1 == null) {
					linkExportString += this.hinges[i].exportString(index);
				} else if (this.hinges[i].body1Index != -1) {
					linkExportString += this.hinges[i].exportString(index);
				}
			}
		}

		for (var i=0; i<this.springs.length; i++) {
			if (this.springs[i].body1 == this) {
				this.springs[i].body1Index = index;
				if (this.springs[i].body2 == null) {
					linkExportString += this.springs[i].exportString(index);
				} else if (this.springs[i].body2Index != -1) {
					linkExportString += this.springs[i].exportString(index);
				}
			} else if (this.springs[i].body2 == this) {
				this.springs[i].body2Index = index;
				if (this.springs[i].body1 == null) {
					linkExportString += this.springs[i].exportString(index);
				} else if (this.springs[i].body1Index != -1) {
					linkExportString += this.springs[i].exportString(index);
				}
			}
		}

		let partOfLevelString = ',0';
		if (this.partOfLevel) {
			partOfLevelString = ',1';
		}

		return shapeExportString + collisionLayerString + partOfLevelString + linkExportString;
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
		this.backgroundColor = 'rgba(40, 40, 40, 1)';

		this.accessToken1 = 'sl.BnFrUs8i4bD-4ZixiSF9Rgg8aDcwhy1zhSrh9qJLb_UAl9rCsEjGjN1BowDQ';
		this.accessToken2 = 'j-3ba1hVZ1PrKsc0jkIjji7oZoowluwiVWb-O8THBm1w8MzbGAeiiyhIGVBL7X1v_ryronkmZYRVDIpa';

		this.home = true;
		this.homeMode = 0;
		this.leaderboardLevel = 0;
		this.level = 0;
		this.numLevels = 10;
		this.username = 'Anonymous';
		this.currLevelStr = '';
		this.playerLevelStr = '';
		this.topScores = {};

		this.startBox = null;
		this.goal = null;
		this.goldenBox = null;
		this.victory = false;

		this.startBoxMinX = 0;
		this.startBoxMaxX = 0;
		this.startBoxMinY = 0;
		this.startBoxMaxY = 0;

		this.victoryScreenDestination = 350;
		this.victoryScreenY = null;
		this.victoryScreenHeight = 400;
		this.victoryScreenWidth = 500;
		this.victoryScreenCurveSize = 16;

		this.ticksPerSecond = 120;
		this.tickStartLevel = 0;
		this.levelTicks = 0;
		this.playerParts = 0;

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
		this.contextMenuHeight = 435;
		this.contextMenuCurveSize = 16;
		this.contextMenuSliderWidth = 0.75 * this.contextMenuWidth;
		this.contextMenuSliderThickness = 3;
		this.contextMenuSliderGrabberWidth = 10;
		this.contextMenuSliderGrabberHeight = 14;
		this.contextMenuSliderGap = 45;
		this.bodySliders = ['angle', 'density', 'restitution', 'friction', 'red', 'green', 'blue', 'alpha'];
		this.bodySlidersMin = [0, 1, 0, 0, 0, 0, 0, 0];
		this.bodySlidersMax = [Math.round(2*Math.PI*100)/100, 5, 1, 1, 255, 255, 255, 1];
		this.rivetSliders = ['size', 'angle', 'red', 'green', 'blue', 'alpha'];
		this.rivetSlidersMin = [5, 0, 0, 0, 0, 0];
		this.rivetSlidersMax = [100, Math.round(2*Math.PI*100)/100, 255, 255, 255, 1];
		this.hingeSliders = ['motorSpeed', 'size', 'red', 'green', 'blue'];
		this.hingeSlidersMin = [-10, 5, 0, 0, 0];
		this.hingeSlidersMax = [10, 100, 255, 255, 255];
		this.springSliders = ['strength', 'damping', 'targetLength', 'width', 'red', 'green', 'blue'];
		this.springSlidersMin = [1, 0, 1, 3, 0, 0, 0];
		this.springSlidersMax = [20, 2, 400, 20, 255, 255, 255];
		this.collisionLayerCheckBoxDistance = this.contextMenuWidth * 0.037;
		this.collisionLayerCheckBoxSize = this.contextMenuWidth/15;

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
		this.levelEditor = false;
	}

	getScores(level, callback) {
		let game = this;
		let xhr = new XMLHttpRequest();
		xhr.open('GET', 'https://content.dropboxapi.com/2/files/download', true);
		xhr.setRequestHeader('Content-Type', 'application/octet-stream');
		xhr.setRequestHeader('Authorization', 'Bearer ' + this.accessToken1 + this.accessToken2);
		xhr.setRequestHeader('Dropbox-API-Arg', '{\"path\": \"/level' + level + '.txt\"}');
		xhr.onreadystatechange = function() {
			if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
				game.topScores[level] = [[], []];
				if (this.responseText.includes(';')) {
					let timeRecords = this.responseText.split(';')[0].split(',');
					let partsRecords = this.responseText.split(';')[1].split(',');
					for (var i=0; i<timeRecords.length; i+=3) {
						game.topScores[level][0][i/3] = [timeRecords[i], parseInt(timeRecords[i+1]), parseInt(timeRecords[i+2])];
					}
					for (var i=0; i<partsRecords.length; i+=3) {
						game.topScores[level][1][i/3] = [partsRecords[i], parseInt(partsRecords[i+1]), parseInt(partsRecords[i+2])];
					}
				}

				callback();
			}
		}
		xhr.send(null);
	}

	resetLevel() {
		this.resetVariables();
		this.parseLevelString(this.currLevelStr + this.playerLevelStr);
	}

	goHome() {
		this.home = true;
		this.homeMode = 0;
		this.resetVariables();
		this.camera.x = 0;
		this.camera.y = 0;
		this.camera.zoom = 1;
		this.currLevelStr = '';
		this.playerLevelStr = '';
	}

	resetVariables() {
		this.mousePosition = new Vector(0, 0);
		this.oldMousePosition = new Vector(0, 0);
		this.contextMenu = null;
		this.dragging = null;
		this.placing = null;
		this.placingDrag = false;
		this.draggingSlider = null;
		this.victory = false;
		this.tickStartLevel = 0;
		this.levelTicks = 0;
		this.pause = true;
		this.canInteract = true;
		this.collisions = [];
		this.bodies = [];
		this.startBox = null;
		this.goal = null;
		this.goldenBox = null;
		this.currid = 0;
		this.playerParts = 0;
	}

	startNextLevel() {
		let game = this;
		this.level++;
		this.loadLevel(this.level + '.lvl', function(levelText) {
			game.resetVariables();
			game.camera.x = 0;
			game.camera.y = 0;
			game.camera.zoom = 1;
			game.currLevelStr = levelText;
			game.playerLevelStr = '';
			game.parseLevelString(levelText);
		});
	}

	calibrateStartBox() {
		if (this.startBox != null) {
			this.startBoxMinX = this.startBox.vertices[0].x;
			this.startBoxMaxX = this.startBox.vertices[0].x;
			this.startBoxMinY = this.startBox.vertices[0].y;
			this.startBoxMaxY = this.startBox.vertices[0].y;

			for (var i=0; i<this.startBox.vertices.length; i++) {
				if (this.startBox.vertices[i].x < this.startBoxMinX) {
					this.startBoxMinX = this.startBox.vertices[i].x;
				}
				if (this.startBox.vertices[i].x > this.startBoxMaxX) {
					this.startBoxMaxX = this.startBox.vertices[i].x;
				}
				if (this.startBox.vertices[i].y < this.startBoxMinY) {
					this.startBoxMinY = this.startBox.vertices[i].y;
				}
				if (this.startBox.vertices[i].y > this.startBoxMinY) {
					this.startBoxMaxY = this.startBox.vertices[i].y;
				}
			}
		}
	}

	parseLevelString(input) {
		let inputArray = input.split(',');

		if (inputArray[0] != 'l' && inputArray[0] != 's' && inputArray[0] != 'h' && inputArray[0] != 'c' && inputArray[0] != 'r') {
			this.backgroundColor = 'rgba(' + inputArray[0] + ', ' + inputArray[1] + ', ' + inputArray[2] + ', ' + inputArray[3] + ')';
			inputArray.splice(0, 4);
		}

		let tempBody;
		let tempLink;
		let links = [];
		let i=0;
		while (i < inputArray.length) {
			switch (inputArray[i]) {
				case 'st':
					this.startBox = new Rectangle([new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])), new Vector(Number(inputArray[i+3]), Number(inputArray[i+4])), 
						new Vector(Number(inputArray[i+5]), Number(inputArray[i+6])), new Vector(Number(inputArray[i+7]), Number(inputArray[i+8]))]);

					this.startBox.collisionLayers = {0: false, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};
					this.startBox.partOfLevel = true;
					this.calibrateStartBox()
					i += 9;
					break;
				case 'g':
					this.goal = new Rectangle([new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])), new Vector(Number(inputArray[i+3]), Number(inputArray[i+4])), 
						new Vector(Number(inputArray[i+5]), Number(inputArray[i+6])), new Vector(Number(inputArray[i+7]), Number(inputArray[i+8]))]);

					this.goal.collisionLayers = {0: false, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};
					this.goal.partOfLevel = true;
					i += 9;
					break;
				case 'b':
					if (this.goldenBox != null) {
						this.removeBody(this.goldenBox);
					}

					this.goldenBox = new Rectangle([new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])), new Vector(Number(inputArray[i+3]), Number(inputArray[i+4])), 
						new Vector(Number(inputArray[i+5]), Number(inputArray[i+6])), new Vector(Number(inputArray[i+7]), Number(inputArray[i+8]))]);

					this.goldenBox.collisionLayers = {0: true, 1: true, 2: true, 3: true, 4: true, 5: true, 6: true, 7: true, 8: true, 9: true};
					this.goldenBox.partOfLevel = true;
					this.goldenBox.goal = true;
					this.addBody(this.goldenBox);
					i += 9;
					break;
				case 'c':
					tempBody = new Circle(new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])), Number(inputArray[i+3]));
					tempBody.angle = Number(inputArray[i+4]);
					tempBody.density = Number(inputArray[i+5]);
					tempBody.restitution = Number(inputArray[i+6]);
					tempBody.friction = Number(inputArray[i+7]);
					tempBody.red = Number(inputArray[i+8]);
					tempBody.green = Number(inputArray[i+9]);
					tempBody.blue = Number(inputArray[i+10]);
					tempBody.alpha = Number(inputArray[i+11]);

					tempBody.collisionLayers = {0: false, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};
					for (var j=0; j<inputArray[i+12].length; j++) {
						tempBody.collisionLayers[Number(inputArray[i+12].charAt(j))] = true;
					}

					tempBody.computeArea();
					tempBody.computeMass();

					tempBody.partOfLevel = (1 == Number(inputArray[i+13]));
					if (!tempBody.partOfLevel) {
						this.playerParts++;
					}

					this.addBody(tempBody);
					i += 14;
					break;
				case 'r':
					tempBody = new Rectangle([new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])), new Vector(Number(inputArray[i+3]), Number(inputArray[i+4])), 
						new Vector(Number(inputArray[i+5]), Number(inputArray[i+6])), new Vector(Number(inputArray[i+7]), Number(inputArray[i+8]))]);

					tempBody.angle = Number(inputArray[i+9]);
					tempBody.density = Number(inputArray[i+10]);
					tempBody.restitution = Number(inputArray[i+11]);
					tempBody.friction = Number(inputArray[i+12]);
					tempBody.red = Number(inputArray[i+13]);
					tempBody.green = Number(inputArray[i+14]);
					tempBody.blue = Number(inputArray[i+15]);
					tempBody.alpha = Number(inputArray[i+16]);

					tempBody.collisionLayers = {0: false, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};
					for (var j=0; j<inputArray[i+17].length; j++) {
						tempBody.collisionLayers[Number(inputArray[i+17].charAt(j))] = true;
					}

					tempBody.computeArea();
					tempBody.computeMass();
					tempBody.computeNormals();

					tempBody.partOfLevel = (1 == Number(inputArray[i+18]));
					if (!tempBody.partOfLevel) {
						this.playerParts++;
					}

					this.addBody(tempBody);
					i += 19;
					break;
				case 'l':
					tempLink = new Rivet(new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])));
					tempLink.body1Index = Number(inputArray[i+3]);
					tempLink.body2Index = Number(inputArray[i+4]);
					tempLink.size = Number(inputArray[i+5]);
					tempLink.angle = Number(inputArray[i+6]);
					tempLink.red = Number(inputArray[i+7]);
					tempLink.green = Number(inputArray[i+8]);
					tempLink.blue = Number(inputArray[i+9]);
					tempLink.alpha = Number(inputArray[i+10]);
					tempLink.partOfLevel = (1 == Number(inputArray[i+11]));
					if (!tempLink.partOfLevel) {
						this.playerParts++;
					}

					links.push(tempLink);
					i += 12;
					break;
				case 's':
					tempLink = new Spring(new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])));
					tempLink.position2.set(Number(inputArray[i+3]), Number(inputArray[i+4]));
					tempLink.body1Index = Number(inputArray[i+5]);
					tempLink.body2Index = Number(inputArray[i+6]);
					tempLink.strength = Number(inputArray[i+7]);
					tempLink.damping = Number(inputArray[i+8]);
					tempLink.width = Number(inputArray[i+9]);
					tempLink.targetLength = Number(inputArray[i+10]);
					tempLink.red = Number(inputArray[i+11]);
					tempLink.green = Number(inputArray[i+12]);
					tempLink.blue = Number(inputArray[i+13]);
					tempLink.partOfLevel = (1 == Number(inputArray[i+14]));
					if (!tempLink.partOfLevel) {
						this.playerParts++;
					}

					links.push(tempLink);
					i += 15;
					break;
				case 'h':
					tempLink = new Hinge(new Vector(Number(inputArray[i+1]), Number(inputArray[i+2])));
					tempLink.body1Index = Number(inputArray[i+3]);
					tempLink.body2Index = Number(inputArray[i+4]);
					tempLink.motorSpeed = Number(inputArray[i+5]);
					tempLink.size = Number(inputArray[i+6]);
					tempLink.red = Number(inputArray[i+7]);
					tempLink.green = Number(inputArray[i+8]);
					tempLink.blue = Number(inputArray[i+9]);
					tempLink.partOfLevel = (1 == Number(inputArray[i+10]));
					if (!tempLink.partOfLevel) {
						this.playerParts++;
					}

					links.push(tempLink);
					i += 11;
					break;
			}
		}

		for (var index=0; index<links.length; index++) {
			this.addLink(links[index]);

			links[index].body1Index = -1;
			links[index].body2Index = -1;
		}
	}

	loadLevel(level, func) {
		return new Promise(function (resolve, reject) {
			let xhr = new XMLHttpRequest();
			xhr.open('GET', 'levels/' + level);
			xhr.onreadystatechange = function() {
				if (this.readyState === XMLHttpRequest.DONE) {
					resolve(this.responseText);
				}
			}
			xhr.send();
		}).then(function(values) {
			if (func) {func(values);}
		}, function() {
			throw ('Error loading level \"levels/' + level + '\"');
		});
	}

	exportPlayer() {
		let exportString = '';
		for (let i=0; i<this.bodies.length; i++) {
			if (!this.bodies[i].partOfLevel) {
				exportString += this.bodies[i].exportString(i);
			}
		}

		if (this.goldenBox != null) {
			exportString += this.goldenBox.getShapeExportString(-3);
		}

		return exportString;
	}

	exportLevel() {
		let exportString = this.backgroundColor.split('(')[1].replace(/ /g, '').replace(')', '');

		if (this.startBox != null) {
			exportString += this.startBox.getShapeExportString(-1);
		}
		if (this.goal != null) {
			exportString += this.goal.getShapeExportString(-2);
		}
		if (this.goldenBox != null) {
			exportString += this.goldenBox.getShapeExportString(-3);
		}

		for (let i=0; i<this.bodies.length; i++) {
			if (this.bodies[i] != this.goldenBox && this.bodies[i].partOfLevel) {
				exportString += this.bodies[i].exportString(i);
			}
		}

		return exportString;
	}

	setupListeners() {
		let gameForListeners = this;
		addMouseDownListener(function(which, eventX, eventY) {
			if (gameForListeners.home) {
				if (gameForListeners.homeMode == 0) {
					switch(which) {
						case 1:
							if (eventX >= gameForListeners.canvas.width/2 - 64 && eventX <= gameForListeners.canvas.width/2 + 64 &&
								eventY >= gameForListeners.canvas.height/5 + 265 && eventY <= gameForListeners.canvas.height/5 + 265 + 30) {

								let userInput = prompt('What is your name?');
								if (userInput && userInput != null && userInput.length > 0 && !userInput.includes(';') && !userInput.includes(',')) {
									gameForListeners.username = userInput.substring(0, 16);
								}
							} else if (eventX >= gameForListeners.canvas.width/2 - 150 && eventX <= gameForListeners.canvas.width/2 + 150 &&
								eventY >= gameForListeners.canvas.height/5 + 385 && eventY <= gameForListeners.canvas.height/5 + 485) {

								gameForListeners.home = false;
								gameForListeners.startNextLevel();
							} else if (eventX >= gameForListeners.canvas.width/2 - 100 && eventX <= gameForListeners.canvas.width/2 + 100 &&
								eventY >= gameForListeners.canvas.height/5 + 510 && eventY <= gameForListeners.canvas.height/5 + 560) {

								gameForListeners.homeMode = 2;
							} else if (eventX >= gameForListeners.canvas.width/2 - 100 && eventX <= gameForListeners.canvas.width/2 + 100 &&
								eventY >= gameForListeners.canvas.height/5 + 585 && eventY <= gameForListeners.canvas.height/5 + 635) {

								gameForListeners.homeMode = 1;
							}
							break;
					}
				} else if (gameForListeners.homeMode == 1) {
					if (eventX >= gameForListeners.canvas.width/2 - 600 && eventX <= gameForListeners.canvas.width/2 - 600 + 200 &&
						eventY >= 25 + 90 * (gameForListeners.numLevels-1) && eventY <= 25 + 90 * (gameForListeners.numLevels-1) + 75) {

						gameForListeners.homeMode = 0;
						return;
					}

					for (var i=0; i<gameForListeners.numLevels; i++) {
						if (eventX >= gameForListeners.canvas.width/2 - 300 && eventX <= gameForListeners.canvas.width/2 + 300 &&
							eventY >= 25 + 90 * i && eventY <= 25 + 90 * i + 75) {

							gameForListeners.home = false;
							gameForListeners.level = i;
							gameForListeners.startNextLevel();
							break;
						}
					}
				} else if (gameForListeners.homeMode == 2) {
					if (eventX >= gameForListeners.canvas.width/2 - 600 && eventX <= gameForListeners.canvas.width/2 - 600 + 200 &&
						eventY >= 25 + 90 * (gameForListeners.numLevels-1) && eventY <= 25 + 90 * (gameForListeners.numLevels-1) + 75) {

						gameForListeners.homeMode = 0;
						return;
					}

					for (var i=0; i<gameForListeners.numLevels; i++) {
						if (eventX >= gameForListeners.canvas.width/2 - 300 && eventX <= gameForListeners.canvas.width/2 + 300 &&
							eventY >= 25 + 90 * i && eventY <= 25 + 90 * i + 75) {

							gameForListeners.homeMode = 3;
							gameForListeners.leaderboardLevel = i+1;
							gameForListeners.getScores(gameForListeners.leaderboardLevel, function() {});
							break;
						}
					}
				} else if (gameForListeners.homeMode == 3) {
					if (eventX >= gameForListeners.canvas.width/2 - 710 && eventX <= gameForListeners.canvas.width/2 - 710 + 200 &&
						eventY >= 25 + 90 * (gameForListeners.numLevels-1) && eventY <= 25 + 90 * (gameForListeners.numLevels-1) + 75) {

						gameForListeners.homeMode = 2;
						return;
					}
				}
			} else {
				let eventPosition = new Vector(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y);
				switch(which) {
					case 1:
						if (gameForListeners.levelEditor) {
							if (eventX >= 0 && eventX <= 200 && eventY >= 0 && eventY <= 100) {
								for (var i=0; i<gameForListeners.bodies.length; i++) {
									gameForListeners.bodies[i].partOfLevel = true;

									for (var j=0; j<gameForListeners.bodies[i].rivets.length; j++) {
										gameForListeners.bodies[i].rivets[j].partOfLevel = true;
									}
									for (var j=0; j<gameForListeners.bodies[i].springs.length; j++) {
										gameForListeners.bodies[i].springs[j].partOfLevel = true;
									}
									for (var j=0; j<gameForListeners.bodies[i].hinges.length; j++) {
										gameForListeners.bodies[i].hinges[j].partOfLevel = true;
									}
								}
								console.log(gameForListeners.exportLevel());
								return;
							} else if (eventX >= 200 && eventX <= 400 && eventY >= 0 && eventY <= 100) {
								let eventXZoomScale = eventX / gameForListeners.camera.zoom;
								let eventYZoomScale = eventY / gameForListeners.camera.zoom;
								gameForListeners.placing = new Rectangle([new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y),
															new Vector(eventXZoomScale + gameForListeners.camera.x + 200, eventYZoomScale + gameForListeners.camera.y),
															new Vector(eventXZoomScale + gameForListeners.camera.x + 200, eventYZoomScale + gameForListeners.camera.y + 100),
															new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y + 100)]);
								gameForListeners.placing.alpha = 0.5;
								gameForListeners.startBox = gameForListeners.placing;
								gameForListeners.startBox.collisionLayers = {0: false, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};
								return;
							} else if (eventX >= 400 && eventX <= 600 && eventY >= 0 && eventY <= 100) {
								let eventXZoomScale = eventX / gameForListeners.camera.zoom;
								let eventYZoomScale = eventY / gameForListeners.camera.zoom;
								gameForListeners.placing = new Rectangle([new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y),
															new Vector(eventXZoomScale + gameForListeners.camera.x + 200, eventYZoomScale + gameForListeners.camera.y),
															new Vector(eventXZoomScale + gameForListeners.camera.x + 200, eventYZoomScale + gameForListeners.camera.y + 100),
															new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y + 100)]);
								gameForListeners.placing.alpha = 0.5;
								gameForListeners.goal = gameForListeners.placing;
								gameForListeners.goal.collisionLayers = {0: false, 1: false, 2: false, 3: false, 4: false, 5: false, 6: false, 7: false, 8: false, 9: false};
								return;
							} else if (eventX >= 600 && eventX <= 800 && eventY >= 0 && eventY <= 100) {
								let eventXZoomScale = eventX / gameForListeners.camera.zoom;
								let eventYZoomScale = eventY / gameForListeners.camera.zoom;
								gameForListeners.placing = new Rectangle([new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y),
															new Vector(eventXZoomScale + gameForListeners.camera.x + 30, eventYZoomScale + gameForListeners.camera.y),
															new Vector(eventXZoomScale + gameForListeners.camera.x + 30, eventYZoomScale + gameForListeners.camera.y + 30),
															new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y + 30)]);
								gameForListeners.placing.alpha = 0.5;
								gameForListeners.goldenBox = gameForListeners.placing;
								gameForListeners.goldenBox.collisionLayers = {0: true, 1: true, 2: true, 3: true, 4: true, 5: true, 6: true, 7: true, 8: true, 9: true};
								gameForListeners.goldenBox.goal = true;
								return;
							}
						}

						if (gameForListeners.victory) {
							if (eventX >= (gameForListeners.canvas.width - gameForListeners.victoryScreenWidth)/2 + gameForListeners.victoryScreenCurveSize + 20 &&
								eventX <= (gameForListeners.canvas.width - gameForListeners.victoryScreenWidth)/2 + gameForListeners.victoryScreenCurveSize + 195 &&
								eventY >= gameForListeners.victoryScreenY + gameForListeners.victoryScreenCurveSize + 300 &&
								eventY <= gameForListeners.victoryScreenY + gameForListeners.victoryScreenCurveSize + 375) {

								gameForListeners.goHome();
							} else if (eventX >= (gameForListeners.canvas.width + gameForListeners.victoryScreenWidth)/2 + gameForListeners.victoryScreenCurveSize - 20 - 175 &&
								eventX <= (gameForListeners.canvas.width + gameForListeners.victoryScreenWidth)/2 + gameForListeners.victoryScreenCurveSize - 20 &&
								eventY >= gameForListeners.victoryScreenY + gameForListeners.victoryScreenCurveSize + 300 &&
								eventY <= gameForListeners.victoryScreenY + gameForListeners.victoryScreenCurveSize + 375) {

								if (gameForListeners.level < gameForListeners.numLevels) {
									gameForListeners.startNextLevel();
								}
							}

							return;
						}

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
								if (eventX >= gameForListeners.contextMenu[1].x + gameForListeners.contextMenuCurveSize + 
									(gameForListeners.contextMenuWidth - gameForListeners.contextMenuSliderWidth - gameForListeners.contextMenuSliderGrabberWidth)/2 &&
									eventX <= gameForListeners.contextMenu[1].x + gameForListeners.contextMenuCurveSize + 
									(gameForListeners.contextMenuWidth - gameForListeners.contextMenuSliderWidth + gameForListeners.contextMenuSliderGrabberWidth)/2 + gameForListeners.contextMenuSliderWidth) {

									let sliders = [];
									if (gameForListeners.contextMenu[0] instanceof Body) {
										sliders = gameForListeners.bodySliders;
									} else if (gameForListeners.contextMenu[0] instanceof Rivet) {
										sliders = gameForListeners.rivetSliders;
									} else if (gameForListeners.contextMenu[0] instanceof Hinge) {
										sliders = gameForListeners.hingeSliders;
									} else if (gameForListeners.contextMenu[0] instanceof Spring) {
										sliders = gameForListeners.springSliders;
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
											} else if (gameForListeners.contextMenu[0] instanceof Link) {
												gameForListeners.contextMenu[0].delete();
												gameForListeners.contextMenu = null;
												gameForListeners.playerParts--;
											}
										}
									} else if (gameForListeners.contextMenu[0] instanceof Body) {
										if (eventY >= gameForListeners.contextMenu[1].y - (gameForListeners.contextMenuHeight - gameForListeners.contextMenuCurveSize- 22 -
												gameForListeners.contextMenuSliderGap * gameForListeners.bodySliders.length) - 5 &&
											eventY <= gameForListeners.contextMenu[1].y - (gameForListeners.contextMenuHeight - gameForListeners.contextMenuCurveSize - 22 -
												gameForListeners.contextMenuSliderGap * gameForListeners.bodySliders.length) + gameForListeners.collisionLayerCheckBoxSize - 5) {

											for (var i=0; i<10; i++) {
												if (eventX >= gameForListeners.contextMenu[1].x + gameForListeners.contextMenuCurveSize + i * (gameForListeners.collisionLayerCheckBoxDistance + gameForListeners.collisionLayerCheckBoxSize) &&
													eventX <= gameForListeners.contextMenu[1].x + gameForListeners.contextMenuCurveSize + i * (gameForListeners.collisionLayerCheckBoxDistance + gameForListeners.collisionLayerCheckBoxSize) + gameForListeners.collisionLayerCheckBoxSize) {

													gameForListeners.contextMenu[0].collisionLayers[i] = !gameForListeners.contextMenu[0].collisionLayers[i];
													break;
												}
											}
										}
									}
								}
							}
						}

						if (gameForListeners.canInteract && !inContextMenu) {
							if (gameForListeners.placing != null) {
								if (gameForListeners.startBox != null && !gameForListeners.levelEditor) {
									if (eventPosition.x > gameForListeners.startBoxMinX &&
										eventPosition.x < gameForListeners.startBoxMaxX &&
										eventPosition.y > gameForListeners.startBoxMinY &&
										eventPosition.y < gameForListeners.startBoxMaxY) {

										gameForListeners.placingDrag = true;
									}
								} else {
									gameForListeners.placingDrag = true;
								}
							} else {
								if (eventX < gameForListeners.buildMenuWidth + gameForListeners.buildMenuCurveSize && 
									eventY < gameForListeners.buildMenuY + 2*gameForListeners.buildMenuCurveSize + gameForListeners.buildMenuHeight && eventY > gameForListeners.buildMenuY) {
									if (eventX < (gameForListeners.buildMenuWidth + gameForListeners.buildMenuCurveSize/1.25 + gameForListeners.buildMenuIconSize)/2 && 
										eventX > (gameForListeners.buildMenuWidth + gameForListeners.buildMenuCurveSize/1.25 - gameForListeners.buildMenuIconSize)/2) {
										for (var i=0; i<5; i++) {
											if (eventY < gameForListeners.buildMenuY + gameForListeners.buildMenuCurveSize + i*(gameForListeners.buildMenuIconSize + 1.5*gameForListeners.buildMenuCurveSize) + gameForListeners.buildMenuIconSize && 
												eventY > gameForListeners.buildMenuY + gameForListeners.buildMenuCurveSize + i*(gameForListeners.buildMenuIconSize + 1.5*gameForListeners.buildMenuCurveSize)) {
												switch(i) {
													case 0:
														gameForListeners.placing = new Circle(eventPosition.copy(), 50);
														gameForListeners.placing.alpha = 0.5;
														break;
													case 1:
														let eventXZoomScale = eventX / gameForListeners.camera.zoom;
														let eventYZoomScale = eventY / gameForListeners.camera.zoom;
														gameForListeners.placing = new Rectangle([new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y),
																					new Vector(eventXZoomScale + gameForListeners.camera.x + 200, eventYZoomScale + gameForListeners.camera.y),
																					new Vector(eventXZoomScale + gameForListeners.camera.x + 200, eventYZoomScale + gameForListeners.camera.y + 100),
																					new Vector(eventXZoomScale + gameForListeners.camera.x, eventYZoomScale + gameForListeners.camera.y + 100)]);
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
													case 4:
														gameForListeners.placing = new Spring(eventPosition.copy());
														gameForListeners.placing.alpha = 0.5;
														break;
												}
											}
										}
									}
								} else {
									if (gameForListeners.startBox == null || gameForListeners.levelEditor ||
										(eventPosition.x > gameForListeners.startBoxMinX &&
										eventPosition.x < gameForListeners.startBoxMaxX &&
										eventPosition.y > gameForListeners.startBoxMinY &&
										eventPosition.y < gameForListeners.startBoxMaxY)) {

										for (var i=gameForListeners.bodies.length-1; i>=0; i--) {
											let hitLink = false;
											for (var j=gameForListeners.bodies[i].rivets.length-1; j>=0; j--) {
												if (gameForListeners.bodies[i].rivets[j].isInside(eventPosition)) {
													gameForListeners.mousePosition.set(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y);
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
													gameForListeners.mousePosition.set(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y);
													gameForListeners.dragging = gameForListeners.bodies[i].hinges[j];
													hitLink = true;
													break;
												}
											}

											if (hitLink) {
												break;
											}

											if (gameForListeners.bodies[i].isInside(eventPosition)) {
												gameForListeners.mousePosition.set(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y);
												gameForListeners.dragging = gameForListeners.bodies[i];
												break;
											}
										}
									}
								}
							}
						}
						break;
					case 3:
						gameForListeners.contextMenu = null;
						if (gameForListeners.canInteract && gameForListeners.dragging == null && gameForListeners.placing == null &&
							(gameForListeners.startBox == null || gameForListeners.levelEditor ||
							(eventPosition.x > gameForListeners.startBoxMinX &&
							eventPosition.x < gameForListeners.startBoxMaxX &&
							eventPosition.y > gameForListeners.startBoxMinY &&
							eventPosition.y < gameForListeners.startBoxMaxY))) {

							for (var i=gameForListeners.bodies.length-1; i>=0; i--) {
								let hitLink = false;
								for (var j=gameForListeners.bodies[i].rivets.length-1; j>=0; j--) {
									if (gameForListeners.bodies[i].rivets[j].isInside(eventPosition)) {
										let contextMenuPosition = new Vector(eventX, eventY);
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
										let contextMenuPosition = new Vector(eventX, eventY);
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

								for (var j=gameForListeners.bodies[i].springs.length-1; j>=0; j--) {
									if (gameForListeners.bodies[i].springs[j].isInside(eventPosition)) {
										let contextMenuPosition = new Vector(eventX, eventY);
										let distanceFromTop = contextMenuPosition.y - gameForListeners.contextMenuHeight - 2*gameForListeners.contextMenuCurveSize - 2;
										let distanceFromSide = window.innerWidth - (contextMenuPosition.x + gameForListeners.contextMenuWidth + 2*gameForListeners.contextMenuCurveSize + 2);
										if (distanceFromTop < 0) {
											contextMenuPosition.addVector(new Vector(0, -distanceFromTop));
										}
										if (distanceFromSide < 0) {
											contextMenuPosition.addVector(new Vector(distanceFromSide, 0));
										}

										gameForListeners.contextMenu = [gameForListeners.bodies[i].springs[j], contextMenuPosition];
										hitLink = true;
										break;
									}
								}

								if (hitLink) {
									break;
								}

								if ((gameForListeners.levelEditor || gameForListeners.bodies[i] != gameForListeners.goldenBox) && gameForListeners.bodies[i].isInside(eventPosition)) {
									let contextMenuPosition = new Vector(eventX, eventY);
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
			}
		});

		addMouseUpListener(function(which, eventX, eventY) {
			if (gameForListeners.home) {
			} else {
				switch(which) {
					case 1:
						let eventPosition = new Vector(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y);

						gameForListeners.dragging = null;
						gameForListeners.draggingSlider = null;

						if (gameForListeners.placing != null && gameForListeners.placingDrag) {
							if (gameForListeners.placing instanceof Body && (getDistance(new Vector(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y), gameForListeners.placing.vertices[0]) > 5 ||
								gameForListeners.placing == gameForListeners.goldenBox)) {

								if (gameForListeners.levelEditor && gameForListeners.placing == gameForListeners.startBox) {
									gameForListeners.placing.place();
									gameForListeners.placing = null;
									gameForListeners.placingDrag = false;
									gameForListeners.calibrateStartBox();
								} else if (gameForListeners.levelEditor && gameForListeners.placing == gameForListeners.goal) {
									gameForListeners.placing.place();
									gameForListeners.placing = null;
									gameForListeners.placingDrag = false;
								} else {
									gameForListeners.placing.computeArea();
									gameForListeners.placing.computeMass();
									if (gameForListeners.placing.vertices.length > 1) {
										gameForListeners.placing.computeNormals();
									}

									gameForListeners.addBody(gameForListeners.placing);

									gameForListeners.placing.place();
									gameForListeners.placing = null;
									gameForListeners.placingDrag = false;
									gameForListeners.playerParts++;
								}
							} else if (gameForListeners.placing instanceof Link && (!(gameForListeners.placing instanceof Spring) || getDistance(new Vector(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x,
								eventY / gameForListeners.camera.zoom + gameForListeners.camera.y), gameForListeners.placing.position) > 5)) {

								gameForListeners.addLink(gameForListeners.placing);

								gameForListeners.placing.place();
								gameForListeners.placing = null;
								gameForListeners.placingDrag = false;
								gameForListeners.playerParts++;
							}
						}
						break;
				}
			}
		});

		addMouseMoveListener(function(eventX, eventY) {
			if (gameForListeners.home) {
			} else {
				gameForListeners.mousePosition.set(eventX / gameForListeners.camera.zoom + gameForListeners.camera.x, eventY / gameForListeners.camera.zoom + gameForListeners.camera.y);
			}
		});

		addMouseWheelListener(function(eventDirection) {
			if (gameForListeners.home) {
			} else {
				if (eventDirection == 1 && gameForListeners.camera.zoom > 0.5) {
					gameForListeners.camera.zoom -= 0.1;
					let mousePos = gameForListeners.mousePosition.copy();
					mousePos.subtractVector(new Vector(gameForListeners.camera.x, gameForListeners.camera.y));

					gameForListeners.camera.x -= 0.1 * mousePos.x * (1/gameForListeners.camera.zoom);
					gameForListeners.camera.y -= 0.1 * mousePos.y * (1/gameForListeners.camera.zoom);
				} else if (eventDirection == -1 && gameForListeners.camera.zoom < 2) {
					gameForListeners.camera.zoom += 0.1;
					let mousePos = gameForListeners.mousePosition.copy();
					mousePos.subtractVector(new Vector(gameForListeners.camera.x, gameForListeners.camera.y));

					gameForListeners.camera.x += 0.1 * mousePos.x * (1/gameForListeners.camera.zoom);
					gameForListeners.camera.y += 0.1 * mousePos.y * (1/gameForListeners.camera.zoom);
				}
			}
		});

		addKeyDownListener(function(keyCode) {
			if (gameForListeners.home) {
			} else {
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
		if (!body.partOfLevel) {
			this.playerParts--;
			this.playerParts -= body.rivets.length;
			this.playerParts -= body.springs.length;
			this.playerParts -= body.hinges.length;
		}

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
			this.playerLevelStr = this.exportPlayer();
			this.dragging = null;
			this.contextMenu = null;
			this.draggingSlider = null;
			this.tickStartLevel = this.tickID;
		} else {
			this.resetLevel();
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

		if (this.home) {
		} else {
			if (this.canInteract) {
				if (this.placing != null) {
					if (this.placingDrag) {
						if (this.placing instanceof Body) {
							if (this.placing instanceof Circle) {
								let newRadius = Math.round(getDistance(this.mousePosition, this.placing.vertices[0]));
								if (this.startBox != null && !this.levelEditor) {
									if (this.placing.vertices[0].x - newRadius > this.startBoxMinX &&
										this.placing.vertices[0].x + newRadius < this.startBoxMaxX &&
										this.placing.vertices[0].y - newRadius > this.startBoxMinY &&
										this.placing.vertices[0].y + newRadius < this.startBoxMaxY) {

										this.placing.radius = newRadius;
									}
								} else {
									this.placing.radius = newRadius;
								}
							} else if (this.placing instanceof Rectangle && this.placing != this.goldenBox) {
								let xDiff = this.mousePosition.x - this.placing.vertices[0].x;
								let yDiff = this.mousePosition.y - this.placing.vertices[0].y;

								let oldVertex1 = this.placing.vertices[1].copy();
								let oldVertex2 = this.placing.vertices[2].copy();
								let oldVertex3 = this.placing.vertices[3].copy();

								if (xDiff * yDiff > 0) {
									this.placing.vertices[1].set(this.mousePosition.x, this.placing.vertices[0].y);
									this.placing.vertices[2].set(this.mousePosition.x, this.mousePosition.y);
									this.placing.vertices[3].set(this.placing.vertices[0].x, this.mousePosition.y);
								} else {
									this.placing.vertices[1].set(this.placing.vertices[0].x, this.mousePosition.y);
									this.placing.vertices[2].set(this.mousePosition.x, this.mousePosition.y);
									this.placing.vertices[3].set(this.mousePosition.x, this.placing.vertices[0].y);
								}

								if (this.startBox != null && !this.levelEditor) {
									let illegal = false;
									for (var i=0; i<this.placing.vertices.length; i++) {
										if (this.placing.vertices[i].x < this.startBoxMinX || 
											this.placing.vertices[i].x > this.startBoxMaxX ||
											this.placing.vertices[i].y < this.startBoxMinY ||
											this.placing.vertices[i].y > this.startBoxMaxY) {

											illegal = true;
											break;
										}
									}

									if (illegal) {
										this.placing.vertices[1] = oldVertex1;
										this.placing.vertices[2] = oldVertex2;
										this.placing.vertices[3] = oldVertex3;
									}
								}
							}
						} else {
							if (this.placing instanceof Spring) {
								let oldPosition2 = this.placing.position2.copy();
								this.placing.position2.set(this.mousePosition.x, this.mousePosition.y);

								if (this.startBox != null && !this.levelEditor &&
									(this.placing.position2.x < this.startBoxMinX || 
									this.placing.position2.x > this.startBoxMaxX ||
									this.placing.position2.y < this.startBoxMinY ||
									this.placing.position2.y > this.startBoxMaxY)) {

									this.placing.position2 = oldPosition2;
								}
							}
						}
					} else {
						this.placing.translate(this.mousePosition.differenceVector(this.oldMousePosition));
					}
				} else if (this.dragging != null) {
					let movementVector = this.mousePosition.differenceVector(this.oldMousePosition);
					this.dragging.translate(movementVector);

					if (this.startBox != null && !this.levelEditor) {
						let illegal = false;
						if (this.dragging instanceof Rectangle) {
							for (var i=0; i<this.dragging.vertices.length; i++) {
								if (this.dragging.vertices[i].x < this.startBoxMinX || 
									this.dragging.vertices[i].x > this.startBoxMaxX ||
									this.dragging.vertices[i].y < this.startBoxMinY ||
									this.dragging.vertices[i].y > this.startBoxMaxY) {
									illegal = true;
									break;
								}
							}
						} else if (this.dragging instanceof Circle) {
							if (this.dragging.vertices[0].x - this.dragging.radius < this.startBoxMinX || 
									this.dragging.vertices[0].x + this.dragging.radius > this.startBoxMaxX ||
									this.dragging.vertices[0].y - this.dragging.radius < this.startBoxMinY ||
									this.dragging.vertices[0].y + this.dragging.radius > this.startBoxMaxY) {
								illegal = true;
							}
						} else if (this.dragging instanceof Link) {
							if (this.dragging.position.x < this.startBoxMinX || 
								this.dragging.position.x > this.startBoxMaxX ||
								this.dragging.position.y < this.startBoxMinY ||
								this.dragging.position.y > this.startBoxMaxY) {
								illegal = true;
							}
						}

						if (illegal) {
							movementVector.negate();
							this.dragging.translate(movementVector);
						}
					}

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
									moved = true;
									break;
								}
							}
						}

						if (!moved) {
							for (var i=0; i<this.dragging.springs.length; i++) {
								if ((this.dragging.springs[i].body1 == this.dragging && !this.dragging.isInside(this.dragging.springs[i].position)) || 
									(this.dragging.springs[i].body2 == this.dragging && !this.dragging.isInside(this.dragging.springs[i].position2))) {

									movementVector.negate();
									this.dragging.translate(movementVector);
									moved = true;
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
					} else if (this.contextMenu[0] instanceof Spring) {
						sliders = this.springSliders;
						slidersMin = this.springSlidersMin;
						slidersMax = this.springSlidersMax;
					}

					let sliderStartX = this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2;
					let percentage = ((this.mousePosition.x - this.camera.x) * this.camera.zoom - sliderStartX)/this.contextMenuSliderWidth;
					let actualNumber = Math.min(Math.max(percentage, 0), 1) * (slidersMax[this.draggingSlider] - slidersMin[this.draggingSlider]) + slidersMin[this.draggingSlider];
					if (Math.abs(slidersMax[this.draggingSlider]) <= 10 && sliders[this.draggingSlider] != 'motorSpeed') {
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
						let stop = true;
						for (var k in this.bodies[i].collisionLayers) {
							if (this.bodies[i].collisionLayers[k] && this.bodies[j].collisionLayers[k]) {
								stop = false;
								break;
							}
						}

						if (stop) {
							continue;
						}

						stop = false;
						for (var k=0; k<this.bodies[i].hinges.length; k++) {
							if (this.bodies[i].hinges[k].getPartner(this.bodies[i]) == this.bodies[j]) {
								stop = true;
								break;
							}
						}

						if (stop) {
							continue;
						}

						if ((this.bodies[i].compositeObject == null && this.bodies[j].compositeObject == null) || this.bodies[i].compositeObject != this.bodies[j].compositeObject) {
							let collision = new Collision(this.bodies[i], this.bodies[j]);
							collision.solve();

							if (collision.contacts.length > 0) {
								this.collisions.push(collision);
							}
						}
					}

					this.bodies[i].applyHingeImpulses();
					this.bodies[i].applySpringImpulses();
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

				if (!this.levelEditor && !this.victory && this.goldenBox != null && this.goal != null) {
					for (var i=0; i<this.goldenBox.vertices.length; i++) {
						if (this.goal.isInside(this.goldenBox.vertices[i])) {
							this.victory = true;
							this.levelTicks = this.tickID - this.tickStartLevel;

							let game = this;
							this.getScores(this.level, function() {
								let update = false;
								if (game.topScores[game.level] && game.topScores[game.level][0] && game.topScores[game.level][0].length > 0) {
									for (var j=0; j<game.topScores[game.level][0].length; j++) {
										if (game.username == game.topScores[game.level][0][j][0]) {
											break;
										}

										if (game.levelTicks < game.topScores[game.level][0][j][1]) {
											game.topScores[game.level][0].splice(j, 0, [game.username, game.levelTicks, game.playerParts]);
											update = true;
											break;
										}

										if (j == game.topScores[game.level][0].length - 1 && j < 10) {
											game.topScores[game.level][0].push([game.username, game.levelTicks, game.playerParts]);
											update = true;
											break;
										}
									}
								} else {
									if (game.topScores[game.level] === undefined) {
										game.topScores[game.level] = [[[game.username, game.levelTicks, game.playerParts]], [[game.username, game.levelTicks, game.playerParts]]];
									}

									if (game.topScores[game.level][0] === undefined || game.topScores[game.level][0].length == 0) {
										game.topScores[game.level][0] = [[game.username, game.levelTicks, game.playerParts]];
									}

									update = true;
								}

								if (game.topScores[game.level] && game.topScores[game.level][1] && game.topScores[game.level][1].length > 0) {
									for (var j=0; j<game.topScores[game.level][1].length; j++) {
										if (game.username == game.topScores[game.level][1][j][0]) {
											break;
										}

										if (game.playerParts < game.topScores[game.level][1][j][2]) {
											game.topScores[game.level][1].splice(j, 0, [game.username, game.levelTicks, game.playerParts]);
											update = true;
											break;
										}

										if (j == game.topScores[game.level][1].length - 1 && j < 10) {
											game.topScores[game.level][1].push([game.username, game.levelTicks, game.playerParts]);
											update = true;
											break;
										}
									}
								} else {
									if (game.topScores[game.level] === undefined) {
										game.topScores[game.level] = [[[game.username, game.levelTicks, game.playerParts]], [[game.username, game.levelTicks, game.playerParts]]];
									}

									if (game.topScores[game.level][1] === undefined || game.topScores[game.level][1].length == 0) {
										game.topScores[game.level][1] = [[game.username, game.levelTicks, game.playerParts]];
									}

									update = true;
								}

								if (update) {
									console.log(game.topScoresToText(game.level));
									let xhr = new XMLHttpRequest();
									xhr.open('POST', 'https://content.dropboxapi.com/2/files/upload', true);
									xhr.setRequestHeader('Content-Type', 'application/octet-stream');
									xhr.setRequestHeader('Authorization', 'Bearer ' + game.accessToken1 + game.accessToken2);
									xhr.setRequestHeader('Dropbox-API-Arg', '{\"path\": \"/level' + game.level + '.txt\",\"mode\": \"overwrite\",\"autorename\": true,\"mute\": false}');
									xhr.onreadystatechange = function() {
										if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
											//console.log('New record!');
										}
									}
									xhr.send(game.topScoresToText(game.level));
								}
							});
							break;
						}
					}
				}

				if (this.victory) {
					if (this.victoryScreenY == null) {
						this.victoryScreenY = this.canvas.height + 10;
					}

					if (Math.abs(this.victoryScreenDestination - this.victoryScreenY) > 3) {
						this.victoryScreenY *= 0.95;
					} else {
						this.victoryScreenY = this.victoryScreenDestination;
					}
				}
			}

			this.oldMousePosition.set(this.mousePosition.x, this.mousePosition.y);
		}
	}

	topScoresToText(level) {
		let scores = this.topScores[level];
		let scoreString = '';
		for (var i=0; i<Math.min(10, this.topScores[level][0].length); i++) {
			scoreString += this.topScores[level][0][i][0] + ',' + this.topScores[level][0][i][1] + ',' + this.topScores[level][0][i][2] + ',';
		}
		scoreString = scoreString.slice(0, -1);
		scoreString += ';';
		for (var i=0; i<Math.min(10, this.topScores[level][1].length); i++) {
			scoreString += this.topScores[level][1][i][0] + ',' + this.topScores[level][1][i][1] + ',' + this.topScores[level][1][i][2] + ',';
		}
		scoreString = scoreString.slice(0, -1);

		return scoreString;
	}

	render() {
		if (this.home) {
			this.context.fillStyle = 'rgba(40, 40, 40, 1)';
			this.context.fillRect(0, 0, this.canvas.width, this.canvas.height);

			if (this.homeMode == 0) {
				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.font = '96px cursive';
				this.context.textAlign = 'center';
				this.context.textBaseline = 'alphabetic';
				this.context.fillText('Riveting', this.canvas.width/2, this.canvas.height/5);
				this.context.font = '24px cursive';
				this.context.fillText('A Physics Game', this.canvas.width/2, this.canvas.height/5 + 65);
				this.context.font = '32px cursive';
				this.context.fillText('Your Name:', this.canvas.width/2, this.canvas.height/5 + 200);
				this.context.fillText(this.username, this.canvas.width/2, this.canvas.height/5 + 245);
				this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
				this.context.lineWidth = 3;
				this.context.fillRect(this.canvas.width/2 - 64, this.canvas.height/5 + 265, 128, 30);
				this.context.strokeRect(this.canvas.width/2 - 64, this.canvas.height/5 + 265, 128, 30);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = '18px cursive';
				this.context.fillText('Change Name', this.canvas.width/2, this.canvas.height/5 + 286);
				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(this.canvas.width/2 - 150, this.canvas.height/5 + 385, 300, 100);
				this.context.strokeRect(this.canvas.width/2 - 150, this.canvas.height/5 + 385, 300, 100);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = '38px cursive';
				this.context.fillText('Play', this.canvas.width/2, this.canvas.height/5 + 435);
				this.context.font = '16px cursive';
				this.context.fillText('Level ' + ((this.level + 1) % this.numLevels), this.canvas.width/2, this.canvas.height/5 + 465);
				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(this.canvas.width/2 - 100, this.canvas.height/5 + 510, 200, 50);
				this.context.strokeRect(this.canvas.width/2 - 100, this.canvas.height/5 + 510, 200, 50);
				this.context.fillRect(this.canvas.width/2 - 100, this.canvas.height/5 + 585, 200, 50);
				this.context.strokeRect(this.canvas.width/2 - 100, this.canvas.height/5 + 585, 200, 50);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = '22px cursive';
				this.context.fillText('Leaderboards', this.canvas.width/2, this.canvas.height/5 + 543);
				this.context.fillText('Level Select', this.canvas.width/2, this.canvas.height/5 + 618);
			} else if (this.homeMode == 1 || this.homeMode == 2) {
				this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
				this.context.lineWidth = 3;
				this.context.font = '28px cursive';
				this.context.textAlign = 'left';
				this.context.textBaseline = 'alphabetic';
				for (var i=0; i<this.numLevels; i++) {
					this.context.fillStyle = 'rgba(255, 255, 255, 1)';
					this.context.fillRect(this.canvas.width/2 - 300, 25 + 90 * i, 600, 75);
					this.context.strokeRect(this.canvas.width/2 - 300, 25 + 90 * i, 600, 75);
					this.context.fillStyle = 'rgba(0, 0, 0, 1)';
					this.context.fillText('Level ' + (i + 1), this.canvas.width/2 - 280, 25 + 90 * i + 47);
				}

				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(this.canvas.width/2 - 600, 25 + 90 * (this.numLevels-1), 200, 75);
				this.context.strokeRect(this.canvas.width/2 - 600, 25 + 90 * (this.numLevels-1), 200, 75);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Back', this.canvas.width/2 - 580, 25 + 90 * (this.numLevels-1) + 47);
			} else if (this.homeMode == 3) {
				this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
				this.context.lineWidth = 3;
				this.context.font = '24px cursive';
				this.context.textAlign = 'left';
				this.context.textBaseline = 'alphabetic';
				for (var i=0; i<10; i++) {
					this.context.fillStyle = 'rgba(255, 255, 255, 1)';
					this.context.fillRect(this.canvas.width/2 - 460, 25 + 90 * i, 435, 75);
					this.context.strokeRect(this.canvas.width/2 - 460, 25 + 90 * i, 435, 75);
					this.context.fillRect(this.canvas.width/2 + 25, 25 + 90 * i, 435, 75);
					this.context.strokeRect(this.canvas.width/2 + 25, 25 + 90 * i, 435, 75);
					this.context.fillStyle = 'rgba(0, 0, 0, 1)';

					if (this.topScores[this.leaderboardLevel] && this.topScores[this.leaderboardLevel][0] && this.topScores[this.leaderboardLevel][0][i]) {
						this.context.fillText((i+1) + '. ' + this.topScores[this.leaderboardLevel][0][i][0] + '  ' +
							((this.topScores[this.leaderboardLevel][0][i][1] / this.ticksPerSecond).toFixed(2)) + 's  ' + this.topScores[this.leaderboardLevel][0][i][2] + ' Parts', this.canvas.width/2 - 440, 25 + 90 * i + 47);
					} else {
						this.context.fillText((i+1) + '. ', this.canvas.width/2 - 440, 25 + 90 * i + 47);
					}

					if (this.topScores[this.leaderboardLevel] && this.topScores[this.leaderboardLevel][1] && this.topScores[this.leaderboardLevel][1][i]) {
						this.context.fillText((i+1) + '. ' + this.topScores[this.leaderboardLevel][1][i][0] + '  ' +
							this.topScores[this.leaderboardLevel][1][i][2] + ' Parts  ' + ((this.topScores[this.leaderboardLevel][1][i][1] / this.ticksPerSecond).toFixed(2)) + 's', this.canvas.width/2 + 45, 25 + 90 * i + 47);
					} else {
						this.context.fillText((i+1) + '. ', this.canvas.width/2 + 45, 25 + 90 * i + 47);
					}
				}

				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(this.canvas.width/2 - 710, 25 + 90 * (this.numLevels-1), 200, 75);
				this.context.strokeRect(this.canvas.width/2 - 710, 25 + 90 * (this.numLevels-1), 200, 75);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Back', this.canvas.width/2 - 680, 25 + 90 * (this.numLevels-1) + 47);
			}
		} else {
			for (var i=0; i<this.inputs.length; i++) {
				switch(this.inputs[i]) {
					case 'KeyW':
						this.camera.y += -5;
						if (this.mousePosition.y != null) {
							this.mousePosition.y += -5;
						}
						break;
					case 'KeyA':
						this.camera.x += -5;
						if (this.mousePosition.x != null) {
							this.mousePosition.x += -5;
						}
						break;
					case 'KeyS':
						this.camera.y += 5;
						if (this.mousePosition.y != null) {
							this.mousePosition.y += 5;
						}
						break;
					case 'KeyD':
						this.camera.x += 5;
						if (this.mousePosition.x != null) {
							this.mousePosition.x += 5;
						}
						break;
					case 'ArrowRight':
						this.camera.x += 5;
						if (this.mousePosition.x != null) {
							this.mousePosition.x += 5;
						}
						break;
					case 'ArrowLeft':
						this.camera.x += -5;
						if (this.mousePosition.x != null) {
							this.mousePosition.x += -5;
						}
						break;
					case 'ArrowUp':
						this.camera.y += -5;
						if (this.mousePosition.y != null) {
							this.mousePosition.y += -5;
						}
						break;
					case 'ArrowDown':
						this.camera.y += 5;
						if (this.mousePosition.y != null) {
							this.mousePosition.y += 5;
						}
						break;
				}
			}

			this.context.fillStyle = this.backgroundColor;
			this.context.fillRect(0, 0, this.canvas.width, this.canvas.height);

			this.context.scale(this.camera.zoom, this.camera.zoom);
			if (this.startBox != null) {
				this.context.strokeStyle = 'rgba(20, 100, 20, 1)';
				this.context.lineWidth = 12;
				this.context.beginPath();
				this.context.moveTo(this.startBox.vertices[3].x - this.camera.x, this.startBox.vertices[3].y - this.camera.y);
				for (var i=0; i<this.startBox.vertices.length; i++) {
					this.context.lineTo(this.startBox.vertices[i].x - this.camera.x, this.startBox.vertices[i].y - this.camera.y);
				}
				this.context.lineTo(this.startBox.vertices[0].x - this.camera.x, this.startBox.vertices[0].y - this.camera.y);
				this.context.stroke();
				this.context.closePath();

				this.context.strokeStyle = 'rgba(70, 200, 70, 1)';
				this.context.lineWidth = 6;
				this.context.beginPath();
				this.context.moveTo(this.startBox.vertices[3].x - this.camera.x, this.startBox.vertices[3].y - this.camera.y);
				for (var i=0; i<this.startBox.vertices.length; i++) {
					this.context.lineTo(this.startBox.vertices[i].x - this.camera.x, this.startBox.vertices[i].y - this.camera.y);
				}
				this.context.lineTo(this.startBox.vertices[0].x - this.camera.x, this.startBox.vertices[0].y - this.camera.y);
				this.context.stroke();
				this.context.closePath();
			}

			for (var i=0; i<this.bodies.length; i++) {
				this.bodies[i].render(this.context, this.camera);

				for (var j=0; j<this.bodies[i].rivets.length; j++) {
					this.bodies[i].rivets[j].render(this.context, this.camera);
				}

				for (var j=0; j<this.bodies[i].hinges.length; j++) {
					this.bodies[i].hinges[j].render(this.context, this.camera);
				}

				for (var j=0; j<this.bodies[i].springs.length; j++) {
					this.bodies[i].springs[j].render(this.context, this.camera);
				}
			}

			if (this.goal != null) {
				this.context.strokeStyle = 'rgba(150, 100, 0, 1)';
				this.context.fillStyle = 'rgba(255, 215, 0, 0.25)';
				this.context.lineWidth = 12;
				this.context.beginPath();
				this.context.moveTo(this.goal.vertices[3].x - this.camera.x, this.goal.vertices[3].y - this.camera.y);
				for (var i=0; i<this.goal.vertices.length; i++) {
					this.context.lineTo(this.goal.vertices[i].x - this.camera.x, this.goal.vertices[i].y - this.camera.y);
				}
				this.context.lineTo(this.goal.vertices[0].x - this.camera.x, this.goal.vertices[0].y - this.camera.y);
				this.context.fill();
				this.context.stroke();
				this.context.closePath();

				this.context.strokeStyle = 'rgba(255, 215, 0, 1)';
				this.context.lineWidth = 6;
				this.context.beginPath();
				this.context.moveTo(this.goal.vertices[3].x - this.camera.x, this.goal.vertices[3].y - this.camera.y);
				for (var i=0; i<this.goal.vertices.length; i++) {
					this.context.lineTo(this.goal.vertices[i].x - this.camera.x, this.goal.vertices[i].y - this.camera.y);
				}
				this.context.lineTo(this.goal.vertices[0].x - this.camera.x, this.goal.vertices[0].y - this.camera.y);
				this.context.stroke();
				this.context.closePath();
			}
			this.context.scale(1/this.camera.zoom, 1/this.camera.zoom);

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

			this.context.font = 'bold 16px sans-serif';
			this.context.textAlign = 'center';
			this.context.textBaseline = 'alphabetic';
			this.context.fillStyle = 'rgba(0, 0, 0, 1)';
			this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
			this.context.lineWidth = 1;
			this.context.fillText('Spring', (this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize - 2 + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize));
			this.context.fillStyle = 'rgba(255, 255, 255, 1)';
			this.context.fillRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
			this.context.strokeRect((this.buildMenuWidth + this.buildMenuCurveSize/1.25 - this.buildMenuIconSize)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize), this.buildMenuIconSize, this.buildMenuIconSize);
			this.context.fillStyle = 'rgba(0, 0, 0, 1)';
			this.context.strokeStyle = 'rgba(0, 0, 255, 1)';
			this.context.lineWidth = 3;
			this.context.beginPath();
			this.context.moveTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.15);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 - 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.275);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 + 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.275);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 - 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.425);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 + 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.425);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 - 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.575);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 + 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.575);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 - 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.725);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2 + 0.2*this.buildMenuIconSize, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.725);
			this.context.lineTo((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.85);
			this.context.stroke();
			this.context.closePath();
			this.context.beginPath();
			this.context.arc((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.15,
								this.buildMenuIconSize * 0.075, 0, 2*Math.PI, false);
			this.context.fill();
			this.context.closePath();
			this.context.beginPath();
			this.context.arc((this.buildMenuWidth + this.buildMenuCurveSize/1.25)/2, this.buildMenuY + this.buildMenuCurveSize + buildMenuIcons*(this.buildMenuIconSize + 1.5*this.buildMenuCurveSize) + this.buildMenuIconSize*0.85,
								this.buildMenuIconSize * 0.075, 0, 2*Math.PI, false);
			this.context.fill();
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

					this.context.font = '12px sans-serif';
					this.context.fillStyle = 'rgba(0, 0, 0, 1)';
					this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
					this.context.textAlign = 'center';
					this.context.textBaseline = 'alphabetic';
					this.context.fillText('Collision Layers', this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2,
						this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 10);

					let cosPi4 = Math.cos(Math.PI/4);
					let sinPi4 = Math.sin(Math.PI/4);

					for (var i=0; i<10; i++) {
						this.context.fillStyle = 'rgba(255, 255, 255, 1)';
						this.context.lineWidth = 1;
						this.context.beginPath();
						this.context.moveTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize),
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5);
						this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize,
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5);
						this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize,
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) + this.collisionLayerCheckBoxSize - 5);
						this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize),
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) + this.collisionLayerCheckBoxSize - 5);
						this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize),
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5);
						this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize,
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5);
						this.context.fill();
						this.context.stroke();
						this.context.closePath();

						this.context.fillStyle = 'rgba(0, 0, 0, 1)';
						this.context.fillText(i, this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize/2,
							this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) + this.collisionLayerCheckBoxSize + 7);

						this.context.lineWidth = 4;
						if (this.contextMenu[0].collisionLayers[i]) {
							this.context.beginPath();
							this.context.moveTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize/2 - (this.collisionLayerCheckBoxSize/2)*cosPi4,
								this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5 + this.collisionLayerCheckBoxSize/2 - (this.collisionLayerCheckBoxSize/2)*sinPi4);
							this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize/2 + (this.collisionLayerCheckBoxSize/2)*cosPi4,
								this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5 + this.collisionLayerCheckBoxSize/2 + (this.collisionLayerCheckBoxSize/2)*sinPi4);
							this.context.moveTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize/2 + (this.collisionLayerCheckBoxSize/2)*sinPi4,
								this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5 + this.collisionLayerCheckBoxSize/2 - (this.collisionLayerCheckBoxSize/2)*cosPi4);
							this.context.lineTo(this.contextMenu[1].x + this.contextMenuCurveSize + i * (this.collisionLayerCheckBoxDistance + this.collisionLayerCheckBoxSize) + this.collisionLayerCheckBoxSize/2 - (this.collisionLayerCheckBoxSize/2)*sinPi4,
								this.contextMenu[1].y - (this.contextMenuHeight - this.contextMenuCurveSize - 22 - this.contextMenuSliderGap * sliders.length) - 5 + this.collisionLayerCheckBoxSize/2 + (this.collisionLayerCheckBoxSize/2)*cosPi4);
							this.context.stroke();
							this.context.closePath();
						}
					}
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
				} else if (this.contextMenu[0] instanceof Spring) {
					sliders = this.springSliders;
					slidersMin = this.springSlidersMin;
					slidersMax = this.springSlidersMax;

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
					if (sliders[j] == 'size' || sliders[j] == 'width' || sliders[j] == 'red' || sliders[j] == 'green' || sliders[j] == 'blue' || sliders[j] == 'alpha' || (sliders[j] == 'angle' && this.contextMenu[0] instanceof Rivet)) {
						this.context.fillStyle = 'rgba(150, 150, 150, 1)';
						this.context.font = '9px sans-serif';
						this.context.textAlign = 'right';
						this.context.fillText('Visual Only', this.contextMenu[1].x + this.contextMenuCurveSize + (this.contextMenuWidth - this.contextMenuSliderWidth)/2 - 4, this.contextMenu[1].y - sliderHeight + this.contextMenuSliderThickness);
					}
					this.context.fillStyle = 'rgba(0, 0, 0, 1)';
					this.context.font = '12px sans-serif';
					this.context.textAlign = 'center';
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
				this.context.scale(this.camera.zoom, this.camera.zoom);
				this.placing.render(this.context, this.camera);
				this.context.scale(1/this.camera.zoom, 1/this.camera.zoom);
			}

			let ticksSoFar = 0;
			if (this.tickStartLevel != 0) {
				ticksSoFar = this.tickID - this.tickStartLevel;
			}
			if (this.victory) {
				ticksSoFar = this.levelTicks
			}
			this.context.fillStyle = 'rgba(0, 0, 0, 1)';
			this.context.font = '32px sans-serif';
			this.context.textAlign = 'right';
			this.context.textBaseline = 'alphabetic';
			this.context.fillText('Time: ' + Math.floor(ticksSoFar / this.ticksPerSecond) + 's', this.canvas.width - 32, 44);
			this.context.fillText('Parts: ' + this.playerParts, this.canvas.width - 32, 86);

			this.context.font = '20px sans-serif';
			this.context.textAlign = 'center';
			if (this.pause) {
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.strokeStyle = 'rgba(200, 0, 0, 1)';
				this.context.lineWidth = 15;
				this.context.beginPath();
				this.context.moveTo(this.canvas.width/2 - 10, 20);
				this.context.lineTo(this.canvas.width/2 - 10, 60);
				this.context.moveTo(this.canvas.width/2 + 10, 20);
				this.context.lineTo(this.canvas.width/2 + 10, 60);
				this.context.stroke();
				this.context.closePath();
				this.context.fillText('Paused. Press Space to Play', this.canvas.width/2, 90);
			} else {
				this.context.fillStyle = 'rgba(0, 200, 0, 1)';
				this.context.lineWidth = 15;
				this.context.beginPath();
				this.context.moveTo(this.canvas.width/2 - 20, 60);
				this.context.lineTo(this.canvas.width/2 - 20, 20);
				this.context.lineTo(this.canvas.width/2 + 20, 40);
				this.context.lineTo(this.canvas.width/2 - 20, 60);
				this.context.fill();
				this.context.closePath();
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Playing. Press Space to Reset', this.canvas.width/2, 90);
			}

			if (this.level == 1) {
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Get the yellow box into the yellow goal to win!', this.canvas.width/2, 150);
				this.context.fillText('You may only build in the green limited space', this.canvas.width/2, 180);
				this.context.fillText('Right click on objects you\'ve created to change their properties', this.canvas.width/2, 210);
			}

			if (this.victory && this.victoryScreenY != null) {
				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
				this.context.lineWidth = 3;
				this.context.beginPath();
				this.context.moveTo((this.canvas.width - this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY);
				this.context.lineTo((this.canvas.width + this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY);
				this.context.arc((this.canvas.width + this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY + this.victoryScreenCurveSize, this.victoryScreenCurveSize, 3*Math.PI/2, 2*Math.PI);
				this.context.lineTo((this.canvas.width + this.victoryScreenWidth)/2 + 2*this.victoryScreenCurveSize, this.victoryScreenY + this.victoryScreenHeight + this.victoryScreenCurveSize);
				this.context.arc((this.canvas.width + this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY + this.victoryScreenHeight + this.victoryScreenCurveSize, this.victoryScreenCurveSize, 0, Math.PI/2);
				this.context.lineTo((this.canvas.width - this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY + this.victoryScreenHeight + 2*this.victoryScreenCurveSize);
				this.context.arc((this.canvas.width - this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY + this.victoryScreenHeight + this.victoryScreenCurveSize, this.victoryScreenCurveSize, Math.PI/2, Math.PI);
				this.context.lineTo((this.canvas.width - this.victoryScreenWidth)/2, this.victoryScreenY + this.victoryScreenCurveSize);
				this.context.arc((this.canvas.width - this.victoryScreenWidth)/2 + this.victoryScreenCurveSize, this.victoryScreenY + this.victoryScreenCurveSize, this.victoryScreenCurveSize, Math.PI, 3*Math.PI/2);
				this.context.fill();
				this.context.stroke();
				this.context.closePath();

				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = '64px sans-serif';
				this.context.textAlign = 'center';
				this.context.textBaseline = 'alphabetic';
				this.context.fillText('Victory!', this.canvas.width/2 + 14, this.victoryScreenY + this.victoryScreenCurveSize + 58);

				this.context.font = '24px sans-serif';
				this.context.fillText('Time: ' + (this.levelTicks / this.ticksPerSecond).toFixed(2) + 's', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 136);
				this.context.fillText('Parts: ' + this.playerParts, this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 168);

				if (this.topScores[this.level] && this.topScores[this.level][0] && this.topScores[this.level][1]) {
					let printed = false;
					for (var i=0; i<Math.min(10, this.topScores[this.level][0].length); i++) {
						if (this.username == this.topScores[this.level][0][i][0]) {
							if (i == 0) {
								this.context.fillText('1st Place Time!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 220);
							} else if (i == 1) {
								this.context.fillText('2nd Place Time!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 220);
							} else if (i == 2) {
								this.context.fillText('3rd Place Time!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 220);
							} else {
								this.context.fillText((i + 1).toString() + 'th Place Time!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 220);
							}
							printed = true;
							break;
						}
					}

					if (!printed) {
						this.context.fillText('>10th Place Time', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 220);
					}

					printed = false;
					for (var i=0; i<Math.min(10, this.topScores[this.level][1].length); i++) {
						if (this.username == this.topScores[this.level][1][i][0]) {
							if (i == 0) {
								this.context.fillText('1st Place Parts!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 252);
							} else if (i == 1) {
								this.context.fillText('2nd Place Parts!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 252);
							} else if (i == 2) {
								this.context.fillText('3rd Place Parts!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 252);
							} else {
								this.context.fillText((i + 1).toString() + 'th Place Parts!', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 252);
							}
							printed = true;
							break;
						}
					}

					if (!printed) {
						this.context.fillText('>10th Place Parts', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 252);
					}
				} else {
					this.context.fillText('Fetching Leaderboards', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 220);
					this.context.fillText('Fetching Leaderboards', this.canvas.width/2 + 10, this.victoryScreenY + this.victoryScreenCurveSize + 252);
				}

				this.context.strokeRect((this.canvas.width - this.victoryScreenWidth)/2 + this.victoryScreenCurveSize + 20, this.victoryScreenY + this.victoryScreenCurveSize + 300, 175, 75);
				this.context.strokeRect((this.canvas.width + this.victoryScreenWidth)/2 + this.victoryScreenCurveSize - 20 - 175, this.victoryScreenY + this.victoryScreenCurveSize + 300, 175, 75);
				this.context.fillText('Home', (this.canvas.width - this.victoryScreenWidth)/2 + this.victoryScreenCurveSize + 107, this.victoryScreenY + this.victoryScreenCurveSize + 345);
				if (this.level == this.numLevels) {
					this.context.fillText('No More Levels!', (this.canvas.width + this.victoryScreenWidth)/2 + this.victoryScreenCurveSize - 20 - 87, this.victoryScreenY + this.victoryScreenCurveSize + 345);
				} else {
					this.context.fillText('Next Level', (this.canvas.width + this.victoryScreenWidth)/2 + this.victoryScreenCurveSize - 20 - 87, this.victoryScreenY + this.victoryScreenCurveSize + 345);
				}
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

			if (this.levelEditor) {
				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
				this.context.lineWidth = 3;
				this.context.fillRect(0, 0, 200, 100);
				this.context.strokeRect(0, 0, 200, 100);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.font = '24px sans-serif';
				this.context.textAlign = 'center';
				this.context.textBaseline = 'middle';
				this.context.fillText('Export', 100, 50);

				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(200, 0, 200, 100);
				this.context.strokeRect(200, 0, 200, 100);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Start Box', 300, 50);

				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(400, 0, 200, 100);
				this.context.strokeRect(400, 0, 200, 100);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Goal', 500, 50);

				this.context.fillStyle = 'rgba(255, 255, 255, 1)';
				this.context.fillRect(600, 0, 200, 100);
				this.context.strokeRect(600, 0, 200, 100);
				this.context.fillStyle = 'rgba(0, 0, 0, 1)';
				this.context.fillText('Golden Box', 700, 50);
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