'use strict'



//  (+1)       (+3)
//   \\       //
//    \\     //
//     ========
//    //     \\
//   //       \\
//  (+2)       (+4)


function Ro(H) {
  return 1.2 * Math.exp(-H/7500);
}

function getADX(alpha) {
	return {
		Cxa: Math.abs(0.135 * Math.cos(alpha)) + Math.abs(0.27 * Math.sin(alpha)),
		Cya: 0.02 * (alpha * 57.3)
	}
}

const KTH = 4.5
const KPSI = 4.5
const g = 9.81

function setControls(state, ctrl, prog) {
	const deltaTh = (prog.getTh(state[0])- state[6]) - state[4] * 3
	const deltaPsi = prog.getPsi(state[0]) - state[7]
	const Th = state[2]

	const _rev1 = 0.5 + deltaTh * KTH + deltaPsi * KPSI
	const _rev2 = 0.5 + deltaTh * KTH - deltaPsi * KPSI
	const _rev3 = 0.5 - deltaTh * KTH + deltaPsi * KPSI
	const _rev4 = 0.5 - deltaTh * KTH - deltaPsi * KPSI
	
	ctrl.rev1 = Math.min(1, Math.max(0, _rev1)) * Math.max(0, 1 - 5.85*Th)
	ctrl.rev2 = Math.min(1, Math.max(0, _rev2))* Math.max(0, 1 - 5.85*Th)
	ctrl.rev3 = Math.min(1, Math.max(0, _rev3))* Math.max(0, 1 - 5.85*Th)
	ctrl.rev4 = Math.min(1, Math.max(0, _rev4))* Math.max(0, 1 - 5.85*Th)
	console.log(ctrl)
}

function derivatives(state, ctrl, prms) {
	const V = state[1]
	const Th = state[2]
	const Psi = state[3]
	const wAlpha = state[4]
	const wGamma = state[5]
	const alpha = state[6]
	const gamma = state[7]
	const Y = state[9]
	const R = prms.w * (ctrl.rev1 + ctrl.rev2 + ctrl.rev3 + ctrl.rev4)
	const Q = Ro(Y) * V * V * 0.5
	const { Cxa, Cya } = getADX(alpha)

	const CA = Math.cos(alpha)
	const SA = Math.sin(alpha)
	const CG = Math.cos(gamma)
	const SG = Math.sin(gamma)
	const CTH = Math.cos(Th)
	const STH = Math.sin(Th)
	const CPSI = Math.cos(Psi)
	const SPSI = Math.sin(Psi)
	const XA = Cxa * Q * prms.S
	const YA = Cya * Q * prms.S
	const G = prms.m * g
	const MV = prms.m * V
	
	const AX = (-R * SA - XA - G * STH) / prms.m
 	const dTH = (R * CA * CG + YA * CG - G * CTH) / MV
    const dPSI = (R * CA * SG + YA * SG) / MV

	const dAlpha = prms.w * (ctrl.rev1 + ctrl.rev2 - ctrl.rev3 - ctrl.rev4) * prms.lx/prms.Jz
	const dGamma = prms.w * (ctrl.rev1 + ctrl.rev3 - ctrl.rev2 - ctrl.rev4) * prms.ly/prms.Jx
	
	return [
	  0,
	  AX,
	  dTH,
	  dPSI,
	  dAlpha,
	  dGamma,
	  wAlpha,
	  wGamma,
	  V * CTH * CPSI,
	  V * STH,
	  V * CTH * SPSI,
	  0,
	  0,
	  0,
	  0,
	]
}


function integrate(state, prms, prog, dT, tMax) {
	const dT_05 = dT * 0.5
	const currentControls = {
		rev1: 0.4,
		rev2: 0.4,
		rev3: 0.6,
		rev4: 0.6,
	}
	let t = 0
	let i = 0
	const result = [state.slice()]

	while(t < tMax) {
		const _state = result[i].slice() 
		setControls(_state, currentControls, prog)

		const K1 = derivatives(state, currentControls, prms)
		const state_K1 = arrayCombination(_state, K1, 1, dT)
		const K2 = derivatives(state_K1, currentControls, prms)
		const derivSumm = arraySumm(K1, K2)
		
		const state2 = arrayCombination(_state, derivSumm, 1, dT_05)
		t += dT
		state2[0] = t
		state2[11] = currentControls.rev1
		state2[12] = currentControls.rev2
		state2[13] = currentControls.rev3
		state2[14] = currentControls.rev4
		result.push(state2)
		i++
	}

	return result
}

//----------------------------------------
function arraySumm(U, V) {
	const count = U.length
	let result = U.slice()
	for(let i = 0; i < count; i++) {
		result[i] += V[i]
	}
	
	return result	
}
// ------------------------------------------
function arrayCombination(U, V, kU, kV) {
	const count = U.length
	const result = []
	for (let i = 0; i < count; i++) {
		result.push(U[i] * kU + V[i] * kV)
	}
	
	return result
}

{
    const dT = 0.01
	                  //t V  Th Psi wA wG  AoA G  X  Y  Z
	const testState = [0, 5, 0, 0,  0, 0,  0,  0, 0, 500, 0, 0.5, 0.5, 0.5, 0.5]
	const testParams = {
		m: 2.5,
		S: 0.032,
		Jz: 1.75,
		Jx: 0.5,
		w: 15,
		lx: 0.2,
		ly: 0.15
	}

	const testProg = {
		getTh: function(t) {
			return -(17.5 / 57.3)
		},
		getPsi: function(t) {
			return 0
		}
	}
	const res = integrate(testState, testParams, testProg, dT, 45.5)
	let output = ''
	res.forEach(r => {
		const [t, V,  Th, Psi, wA, wG,  AoA, G,  X,  Y,  Z, r1, r2, r3, r4] = r
		const str = [
			t.toFixed(2),
			V.toFixed(1),
			(Th * 57.3).toFixed(2),
			(Psi * 57.3).toFixed(2),
			(wA * 57.3).toFixed(2),
			(wG * 57.3).toFixed(2),
			(AoA * 57.3).toFixed(2),
			(G * 57.3).toFixed(2),
			X.toFixed(1),
			Y.toFixed(1),
			Z.toFixed(1),
			r1.toFixed(2),
			r2.toFixed(2),
			r3.toFixed(2),
			r4.toFixed(2),
		].join('\t')

		output += (str + '\n')
 	})

	console.log(output)
}