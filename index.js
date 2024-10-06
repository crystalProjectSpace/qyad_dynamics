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
		Cxa: Math.abs(0.135 * Math.cos(alpha)) + Math.abs(0.27 * Math.sin(alpha))
		Cya: 0.02 * (alpha * 57.3)
	}
}

function controls(t, state) {
	const deltaTh = thProg - state[2]
	const deltaPsi = psiProg - state[3]
	
}

function derivatives(state, ctrl, prms) {
	const V = state[1]
	const Th = state[2]
	const Psi = state[3]
	const dAlpha = state[4]
	const dGamma = state[5]
	const alpha = state[6]
	const gamma = state[7]
	const Y = state[9]
	const R = w * (ctrl.rev1 + ctrl.rev2 + ctrl.rev3 + ctrl.rev4)
	const Q = Ro(H) * V * V * 0.5
	const { Cxa, Cya } = getADX(alpha)
	const CA = Math.cos(alpha)
	const SA = Math.sin(alpha)
	const CG = Math.cos(gamma)
	const SG = Math.sin(gamma)
	const CTH = Math.cos(Th)
	const STH = Math.sin(Th)
	const CPSI = Math.cos(Psi)
	const SPSI = Math.sin(Psi)
	const dir = alpha > 0 ? -1 : 1
	const XA = Cxa * Q * prm.S
	const YA = Cya * Q * prm.S
	const G = prm.m * g
	const MV = prm.m * V
	
	const AX = (R * CA  * dir - XA - G * STH) / prm.m
 	const dTH = (R * SA * CG + YA * CG - G * CTH) / MV
    const dPSI = (R * SA * SG + YA * SG) / MV
	
	const dAlpha = w * (ctrl.rev1 + ctrl.rev2 - ctrl.rev3 - ctrl.rev4)/prm.Jz
	const dGamma = w * (ctrl.rev1 + ctrl.rev3 - ctrl.rev2 - ctrl.rev4)/prm.Jx
	
	return [
	  AX,
	  dTH,
	  dPSI,
	  dAlpha,
	  dGamma,
	  wZ,
	  wX,
	  V * CTH * CPSI
	  V * STH
	  V * CTH * SPSI
	]
}

function integrate(state, prms) {
	const dT_05 = dT * 0.5
	
	while(t < tMax) {
		const K1 = derivatives(state, ctrl, prms)
		const state_K1 = arrayCombination(state, K1, 1, dT)
		const K2 = derivatives(state_1, ctrl, prms)
		const derivSumm = arraySumm(K1, K2)
		const state2 = arrayCombination(state, derivSumm, 1, dT_05)
	}
}

const testParams = {
	mass: 2.5
	S: 0.032
	
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