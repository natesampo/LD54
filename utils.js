function toDegrees(radians) {
	return radians * (180 / Math.PI);
}

function toRadians(degrees) {
	return degrees * (Math.PI / 180);
}

function contains(list, item) {
	for (var i=0; i<list.length; i++) {
		if (list[i] == item) {
			return true;
		}
	}

	return false;
}

function remove(list, item) {
	for (var i=0; i<list.length; i++) {
		if (list[i] == item) {
			list.splice(i, 1);
			break;
		}
	}
}

function getIndex(list, item) {
	for (var i=0; i<list.length; i++) {
		if (list[i] == item) {
			return i;
		}
	}

	return -1;
}

function copyArray(list) {
	let copy = [];
	for (var i=0; i<list.length; i++) {
		copy.push(list[i]);
	}

	return copy;
}

function getDistance(vector1, vector2) {
	let a = vector2.x - vector1.x;
	let b = vector2.y - vector1.y;
	return Math.sqrt(a * a + b * b);
}

function pointToLineDistance(point, vector1, vector2) {
	let a = vector1.y - vector2.y;
	let b = vector2.x - vector1.x;
	let c = (vector1.x - vector2.x) * vector1.y + (vector2.y - vector1.y) * vector1.x;

	return Math.abs(a * point.x + b * point.y + c) / Math.sqrt(a * a + b * b);
}

function areaOfTriangle(vector1, vector2, vector3) {
	return Math.abs((vector2.x * vector1.y - vector1.x * vector2.y) +
					(vector3.x * vector2.y - vector2.x * vector3.y) +
					(vector1.x * vector3.y - vector3.x * vector1.y)) / 2;
}

function biasedGreaterThan(a, b) {
	let biasRelative = 0.95;
	let biasAbsolute = 0.01;

	return a >= b * biasRelative + a * biasAbsolute;
}

function clip(normal, c, face) {
	let face1 = face[0].copy();
	let face2 = face[1].copy();
	let out = [face1, face2];
	let sp = 0;

	let d1 = normal.dotProduct(face1) - c;
	let d2 = normal.dotProduct(face2) - c;

	if (d1 <= 0) {
		out[sp] = face1;
		sp++;
	}
	if (d2 <= 0) {
		out[sp] = face2;
		sp++;
	}
	if (d1 * d2 < 0) {
		let alpha = d1 / (d1 - d2);
		let newFace2 = face2.copy();
		newFace2.subtractVector(face1);
		newFace2.multiplyScalar(alpha);
		newFace2.addVector(face1)
		out[sp] = newFace2;
		sp++;
	}

	return [sp, out[0], out[1]];
}

function scalarVectorCrossProduct(scalar, vector) {
	return new Vector(-scalar * vector.y, scalar * vector.x);
}

function getAngle(vector1, vector2) {
	let dy = vector2.y - vector1.y;
	let dx = vector2.x - vector1.x;
	let theta = Math.atan2(dy, dx);
	return theta;
}

// takes a point of form vector and an array of connected points of form [vector, vector...]
function pointInPolygon(point, polygon) {
	var x = point.x, y = point.y;

	var inside = false;
	for (var i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
		var xi = polygon[i].x, yi = polygon[i].y;
		var xj = polygon[j].x, yj = polygon[j].y;

		var intersect = ((yi > y) != (yj > y))
		&& (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
		if (intersect) inside = !inside;
	}

	return inside;
};