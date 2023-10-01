class Collision {
	constructor(body1, body2) {
		this.body1 = body1;
		this.body2 = body2;

		this.restitution = Math.min(body1.restitution, body2.restitution);
		this.friction = Math.sqrt(body1.friction * body2.friction);
		this.contacts = [];
		this.penetration = 0;
		this.normal = new Vector(0, 0);
	}

	solve() {
		if (this.body1 instanceof Circle) {
			if (this.body2 instanceof Circle) {
				this.circleToCircle();
			} else if (this.body2 instanceof Rectangle) {
				this.circleToRectangle();
			}
		} else if (this.body1 instanceof Rectangle) {
			if (this.body2 instanceof Circle) {
				let tempBody = this.body1;
				this.body1 = this.body2;
				this.body2 = tempBody;
				this.circleToRectangle();
			} else if (this.body2 instanceof Rectangle) {
				this.rectangleToRectangle();
			}
		}
	}

	applyImpulse() {
		if (this.body1.getMass() + this.body2.getMass() == 0) {
			this.infiniteMassCorrection();
			return;
		}

		for (var i=0; i<this.contacts.length; i++) {
			let distance1 = this.contacts[i].copy();
			distance1.subtractVector(this.body1.findCenterOfMass());

			let distance2 = this.contacts[i].copy();
			distance2.subtractVector(this.body2.findCenterOfMass());

			let relativeVelocity = this.body2.velocity.copy();
			relativeVelocity.addVector(scalarVectorCrossProduct(this.body2.angularVelocity, distance2));
			relativeVelocity.subtractVector(this.body1.velocity);
			relativeVelocity.subtractVector(scalarVectorCrossProduct(this.body1.angularVelocity, distance1));

			let contactVelocity = relativeVelocity.dotProduct(this.normal);
			if (contactVelocity > 0) {
				return;
			}

			let distance1CrossNormal = distance1.vectorCrossProduct(this.normal);
			let distance2CrossNormal = distance2.vectorCrossProduct(this.normal);
			let invMassSum = this.body1.getInvMass() + this.body2.getInvMass() + distance1CrossNormal * distance1CrossNormal * this.body1.getInvMomentOfInertia() +
								distance2CrossNormal * distance2CrossNormal * this.body2.getInvMomentOfInertia();

			let j = -(1 + this.restitution) * contactVelocity;
			j /= invMassSum;
			j /= this.contacts.length;

			let impulse = this.normal.copy();
			impulse.multiplyScalar(j);
			this.body2.applyImpulse(impulse, distance2);
			impulse.negate();
			this.body1.applyImpulse(impulse, distance1);

			relativeVelocity = this.body2.velocity.copy();
			relativeVelocity.addVector(scalarVectorCrossProduct(this.body2.angularVelocity, distance2));
			relativeVelocity.subtractVector(this.body1.velocity);
			relativeVelocity.subtractVector(scalarVectorCrossProduct(this.body1.angularVelocity, distance1));

			let tangent = this.normal.copy();
			tangent.multiplyScalar(relativeVelocity.dotProduct(this.normal));
			tangent.negate();
			tangent.addVector(relativeVelocity);
			tangent.normalize();

			let jTangent = -relativeVelocity.dotProduct(tangent);
			jTangent /= invMassSum;
			jTangent /= this.contacts.length;

			if (jTangent < 0.001 && jTangent > 0.001) {
				return;
			}

			let tangentImpulse = tangent.copy();
			tangentImpulse.multiplyScalar(-j);
			tangentImpulse.multiplyScalar(this.friction);

			this.body2.applyImpulse(tangentImpulse, distance2);
			tangentImpulse.negate();
			this.body1.applyImpulse(tangentImpulse, distance1);
		}
	}

	circleToCircle() {
		this.contacts = [];

		let v = this.body2.vertices[0].copy();
		v.subtractVector(this.body1.vertices[0]);

		let distanceSquared = v.x * v.x + v.y * v.y;
		let radii = this.body1.radius + this.body2.radius;

		if (distanceSquared > radii * radii) {
			return;
		}

		let distance = Math.sqrt(distanceSquared);
		if (distance == 0) {
			this.penetration = this.body1.radius;
			this.normal = new Vector(1, 0);
			this.contacts.push(this.body1.vertices[0]);
		} else {
			this.penetration = radii - distance;
			v.divideScalar(distance);
			this.normal = v;
			let contactPoint = v.copy();
			contactPoint.multiplyScalar(this.body1.radius);
			contactPoint.addVector(this.body1.vertices[0]);
			this.contacts.push(contactPoint);
		}
	}

	circleToRectangle() {
		this.contacts = [];

		let center = this.body1.vertices[0];

		let maxSeparation = null;
		let faceIndex = null;
		for (var i=0; i<this.body2.vertices.length; i++) {
			let tempCenter = center.copy();
			tempCenter.subtractVector(this.body2.vertices[i]);
			let separation = this.body2.normals[i].dotProduct(tempCenter);

			if (separation > this.body1.radius) {
				return;
			}

			if (maxSeparation == null || separation > maxSeparation) {
				maxSeparation = separation;
				faceIndex = i;
			}
		}

		if (maxSeparation < 0.0001) {
			let faceNormal = this.body2.normals[faceIndex].copy();
			faceNormal.negate();
			this.normal = faceNormal;

			let normalCopy = faceNormal.copy();
			normalCopy.multiplyScalar(this.body1.radius);
			normalCopy.addVector(center);
			this.contacts.push(normalCopy);
			this.penetration = this.body1.radius;

			return;
		}

		this.penetration = this.body1.radius - maxSeparation;

		let faceIndexPlus = (faceIndex+1) % this.body2.vertices.length;
		let centerCopy = center.copy();
		centerCopy.subtractVector(this.body2.vertices[faceIndex]);
		let faceCopy = this.body2.vertices[faceIndexPlus].copy();
		faceCopy.subtractVector(this.body2.vertices[faceIndex]);
		let dot1 = centerCopy.dotProduct(faceCopy);

		centerCopy = center.copy();
		centerCopy.subtractVector(this.body2.vertices[faceIndexPlus]);
		faceCopy = this.body2.vertices[faceIndex].copy();
		faceCopy.subtractVector(this.body2.vertices[faceIndexPlus]);
		let dot2 = centerCopy.dotProduct(faceCopy);

		let diffX = 0;
		let diffY = 0;
		let normal = null;

		if (dot1 <= 0) {
			diffX = center.x - this.body2.vertices[faceIndex].x;
			diffY = center.y - this.body2.vertices[faceIndex].y;
			if (diffX * diffX + diffY * diffY > this.body1.radius * this.body1.radius) {
				return;
			}

			normal = this.body2.vertices[faceIndex].copy();
			normal.subtractVector(center);
			normal.normalize();
			this.normal = normal;
			this.contacts.push(this.body2.vertices[faceIndex]);
		} else if (dot2 <= 0) {
			diffX = center.x - this.body2.vertices[faceIndexPlus].x;
			diffY = center.y - this.body2.vertices[faceIndexPlus].y;
			if (diffX * diffX + diffY * diffY > this.body1.radius * this.body1.radius) {
				return;
			}

			normal = this.body2.vertices[faceIndexPlus].copy();
			normal.subtractVector(center);
			normal.normalize();
			this.normal = normal;
			this.contacts.push(this.body2.vertices[faceIndexPlus]);
		} else {
			centerCopy = center.copy();
			centerCopy.subtractVector(this.body2.vertices[faceIndex]);
			normal = this.body2.normals[faceIndex].copy();
			if (centerCopy.dotProduct(normal) > this.body1.radius) {
				return;
			}

			normal.negate();
			this.normal = normal;
			let newContact = normal.copy();
			newContact.multiplyScalar(this.body1.radius);
			newContact.addVector(center);
			this.contacts.push(newContact);
		}
	}

	rectangleToRectangle() {
		this.contacts = [];

		let penetrationA = this.body1.getLeastPenetration(this.body2);
		if (penetrationA[0] >= 0) {
			return;
		}

		let penetrationB = this.body2.getLeastPenetration(this.body1);
		if (penetrationB[0] >= 0) {
			return;
		}

		let referenceBody = null;
		let incidentBody = null;
		let referenceFace = null;
		let flip = false;
		if (biasedGreaterThan(penetrationA[0], penetrationB[0])) {
			referenceBody = this.body1;
			incidentBody = this.body2;
			referenceFace = penetrationA[1];
			flip = false;
		} else {
			referenceBody = this.body2;
			incidentBody = this.body1;
			referenceFace = penetrationB[1];
			flip = true;
		}

		let v1 = referenceBody.vertices[referenceFace];
		let v2 = referenceBody.vertices[(referenceFace+1) % referenceBody.vertices.length];
		let sidePlaneNormal = v2.copy();
		sidePlaneNormal.subtractVector(v1);
		sidePlaneNormal.normalize();

		let referenceFaceNormal = referenceBody.normals[referenceFace].copy();
		let refC = referenceFaceNormal.dotProduct(v1);
		let negativeSide = -sidePlaneNormal.dotProduct(v1);
		let positiveSide = sidePlaneNormal.dotProduct(v2);

		let incidentIndex = referenceBody.findIncidentFace(incidentBody, referenceFace);
		let incidentFace1 = incidentBody.vertices[incidentIndex];
		let incidentFace2 = incidentBody.vertices[(incidentIndex+1) % incidentBody.vertices.length];
		let clipOutput = clip(sidePlaneNormal, positiveSide, [incidentFace1, incidentFace2]);
		if (clipOutput[0] < 2) {
			return;
		}
		incidentFace1 = clipOutput[1];
		incidentFace2 = clipOutput[2];

		sidePlaneNormal.negate();
		clipOutput = clip(sidePlaneNormal, negativeSide, [incidentFace1, incidentFace2]);
		if (clipOutput < 2) {
			return;
		}
		incidentFace1 = clipOutput[1];
		incidentFace2 = clipOutput[2];

		let clippedPoints = 0;
		let separation = referenceFaceNormal.dotProduct(incidentFace1) - refC;
		if (separation <= 0) {
			this.contacts.push(incidentFace1);
			this.penetration = -separation;
			clippedPoints++;
		} else {
			this.penetration = 0;
		}

		separation = referenceFaceNormal.dotProduct(incidentFace2) - refC;
		if (separation <= 0) {
			this.contacts.push(incidentFace2);
			this.penetration -= separation;
			clippedPoints++;

			this.penetration /= clippedPoints;
		}

		if (flip) {
			referenceFaceNormal.negate();
		}
		this.normal = referenceFaceNormal;
	}

	positionalCorrection() {
		if (this.body1.getMass() + this.body2.getMass() > 0) {
			let slack = 0.01;
			let percent = 0.5;

			let correction = (Math.max(this.penetration - slack, 0) / (this.body1.getInvMass() + this.body2.getInvMass())) * percent;
			let correctionVector = this.normal.copy();
			correctionVector.multiplyScalar(correction);
			let correctionVector2 = correctionVector.copy();
			correctionVector.multiplyScalar(this.body1.getInvMass());
			correctionVector.negate();
			correctionVector2.multiplyScalar(this.body2.getInvMass());

			this.body1.translateComposite(correctionVector);
			this.body2.translateComposite(correctionVector2);
		}
	}

	infiniteMassCorrection() {
		this.body1.velocity.set(0, 0);
		this.body2.velocity.set(0, 0);

		this.body1.angularVelocity = 0;
		this.body2.angularVelocity = 0;
	}
}