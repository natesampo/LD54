function addInputs(inputs) {
	document.addEventListener('keydown', function(event) {
		inputs[event.code] = event.key;
	});

	document.addEventListener('keyup', function(event) {
		delete inputs[event.code];
	});

	document.addEventListener('mousedown', function(event) {
		inputs['mouse' + event.which] = true;
	});

	document.addEventListener('mouseup', function(event) {
		delete inputs['mouse' + event.which];
	});
}

function addKeyDownListener(func) {
	document.addEventListener('keydown', function(event) {
		func(event.code);
	});
}

function addKeyUpListener(func) {
	document.addEventListener('keyup', function(event) {
		func(event.code);
	});
}

function addMouseDownListener(func) {
	document.addEventListener('mousedown', function(event) {
		func(event.which, event.clientX, event.clientY);
	});
}

function addMouseUpListener(func) {
	document.addEventListener('mouseup', function(event) {
		func(event.which, event.clientX, event.clientY);
	});
}

function addMouseMoveListener(func) {
	document.addEventListener('mousemove', function(event) {
		func(event.clientX, event.clientY);
	});
}

function addMouseWheelListener(func) {
	document.addEventListener('wheel', function(event) {
		func(Math.sign(event.deltaY));
	});
}

function preventContextMenu() {
	document.addEventListener('contextmenu', function(event) {
		event.preventDefault();
	});
}