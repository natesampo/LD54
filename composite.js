class CompositeObject {
	constructor(objects, links) {
		this.objects = [];
		this.links = [];
		this.connectivity = {};
		this.numConnections = {};

		this.totalMass = 0;
		this.masslessObjects = 0;
		this.masslessLinks = 0;
		this.momentOfInertia = 0;
		this.invMass = 0;
		this.invMomentOfInertia = 0;

		for (var i=0; i<objects.length; i++) {
			this.addObject(objects[i]);
		}

		for (var i=0; i<links.length; i++) {
			this.addLink(links[i]);
		}
	}

	applyImpulse(impulse, contactVector) {
		for (var i=0; i<this.objects.length; i++) {
			let tempImpulse = impulse.copy();
			tempImpulse.multiplyScalar(this.getInvMass());
			this.objects[i].velocity.addVector(tempImpulse);

			let deltaAngularVelocity = this.objects[i].getInvMomentOfInertia() * contactVector.vectorCrossProduct(impulse);
			this.objects[i].angularVelocity += deltaAngularVelocity;
		}
	}

	applyTorque(torque) {
		for (var i=0; i<this.objects.length; i++) {
			this.objects[i].angularVelocity += torque;
		}
	}

	translate(vector) {
		for (var i=0; i<this.objects.length; i++) {
			this.objects[i].translate(vector);
		}
	}

	getMass() {
		if (this.masslessObjects > 0 || this.masslessLinks > 0) {
			return 0;
		}

		return this.totalMass;
	}

	getInvMass() {
		if (this.masslessObjects > 0 || this.masslessLinks > 0 || this.totalMass == 0) {
			return 0;
		}

		return this.invMass;	
	}

	getInvMomentOfInertia() {
		if (this.masslessObjects > 0 || this.masslessLinks > 0 || this.totalMass == 0) {
			return 0;
		}

		return this.invMomentOfInertia;
	}

	calculateMass() {
		this.masslessObjects = 0;
		this.totalMass = 0;
		this.momentOfInertia = 0;
		this.invMass = 0;
		this.invMomentOfInertia = 0;
		for (var i=0; i<this.objects.length; i++) {
			this.totalMass += this.objects[i].mass;
			this.momentOfInertia += this.objects[i].momentOfInertia;

			if (this.objects[i].masslessObjects == 0) {
				this.masslessObjects++;
			}
		}

		if (this.totalMass > 0) {
			this.invMass = 1/this.totalMass;
			this.invMomentOfInertia = 1/this.momentOfInertia;
		}
	}

	findCenterOfMass() {
		let centerOfMass = new Vector(0, 0);
		for (var i=0; i<this.objects.length; i++) {
			let objectCOM = this.objects[i].findCenterOfMassForComposite();
			objectCOM.multiplyScalar(this.objects[i].mass/this.totalMass);
			centerOfMass.addVector(objectCOM);
		}

		return centerOfMass;
	}

	addObject(object) {
		if (object != null) {
			object.compositeObject = this;
			if (!contains(this.objects, object)) {
				this.objects.push(object);

				this.totalMass += object.mass;
				if (this.totalMass > 0) {
					this.invMass = 1/this.totalMass;
				} else {
					this.invMass = 0;
				}
				this.momentOfInertia += object.momentOfInertia;
				this.invMomentOfInertia = 1/this.momentOfInertia;
				if (object.mass == 0) {
					this.masslessObjects++;
				}
			}
		}
	}

	addLink(link) {
		link.compositeObject = this;
		if (!contains(this.links, link)) {
			this.links.push(link);

			let index1 = null;
			let index2 = null;
			if (link.body1 != null) {
				index1 = link.body1.id;
			}
			if (link.body2 != null) {
				index2 = link.body2.id;
			}

			if (!this.connectivity[index1]) {
				this.connectivity[index1] = [];
				this.numConnections[index1] = {};
			}
			if (!this.connectivity[index2]) {
				this.connectivity[index2] = [];
				this.numConnections[index2] = {};
			}

			if (!this.numConnections[index1][index2]) {
				this.numConnections[index1][index2] = 0;
			}
			if (!this.numConnections[index2][index1]) {
				this.numConnections[index2][index1] = 0;
			}

			this.connectivity[index1].push(index2);
			this.connectivity[index2].push(index1);
			this.numConnections[index1][index2]++;
			this.numConnections[index2][index1]++;

			if (index1 == null || index2 == null) {
				this.masslessLinks++;
			}
		}
	}

	removeLink(link) {
		if (contains(this.links, link)) {
			link.compositeObject = null;
			remove(this.links, link);

			if (this.links.length == 0) {
				for (var i=0; i<this.objects.length; i++) {
					this.objects[i].compositeObject = null;
					this.objects[i].rivets = [];
				}

				this.objects = [];
				this.links = [];
				this.connectivity = {};
				this.numConnections = {};

				return;
			}

			let index1 = null;
			let index2 = null;
			if (link.body1 != null) {
				index1 = link.body1.id;
			}
			if (link.body2 != null) {
				index2 = link.body2.id;
			}

			if (index1 == null || index2 == null) {
				this.masslessLinks--;
			}

			this.numConnections[index1][index2]--;
			this.numConnections[index2][index1]--;

			if (this.numConnections[index1][index2] <= 0) {
				delete this.numConnections[index1][index2];
				delete this.numConnections[index2][index1];
				remove(this.connectivity[index1], index2);
				remove(this.connectivity[index2], index1);

				let visited = {};
				this.depthFirstSearch(this.objects[0].id, visited);

				if (Object.keys(visited).length == 1) {
					if (this.objects.length < 3) {
						for (var i=0; i<this.objects.length; i++) {
							this.objects[i].compositeObject = null;
							this.objects[i].rivets = [];
						}

						this.objects = [];
						this.links = [];
						this.connectivity = {};
						this.numConnections = {};

						return;
					} else {
						this.objects[0].compositeObject = null;
						this.objects[0].rivets = [];

						remove(this.objects, this.objects[0]);
					}
				} else {
					let disconnectedNodes = [];
					let disconnectedLinks = [];
					for (var i=this.objects.length-1; i>=0; i--) {
						if (!visited[this.objects[i].id]) {
							for (var j=0; j<this.objects[i].rivets.length; j++) {
								disconnectedLinks.push(this.objects[i].rivets[j]);
								remove(this.links, this.objects[i].rivets[j]);
							}

							disconnectedNodes.push(this.objects.splice(i, 1)[0]);
						}
					}

					if (disconnectedNodes.length == 1) {
						disconnectedNodes[0].compositeObject = null;
						disconnectedNodes[0].rivets = [];
					} else {
						new CompositeObject(disconnectedNodes, disconnectedLinks);
					}
				}
			}
		}

		this.calculateMass();
	}

	merge(compositeObject) {
		for (var i=0; i<compositeObject.objects.length; i++) {
			this.addObject(compositeObject.objects[i]);
			compositeObject.objects[i].compositeObject = this;
		}
		compositeObject.objects = [];

		for (var i=0; i<compositeObject.links.length; i++) {
			this.addLink(compositeObject.links[i]);
			compositeObject.links[i].compositeObject = this;
		}
		compositeObject.links = [];
	}

	depthFirstSearch(id, visited) {
		visited[id] = true;
		for (var i=0; i<this.connectivity[id].length; i++) {
			if (!visited[this.connectivity[id][i]]) {
				this.depthFirstSearch(this.connectivity[id][i], visited);
			}
		}
	}
}