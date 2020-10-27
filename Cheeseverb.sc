// Mono reverb
Cheeseverb1{
	*ar{ |in, fbamount=0.1, highcut=5300, lowcut=100, verbtime=5.0, modfreq=0.1, damp=0.2, modamount=0.15|
		var fb, sig;
		var source = in;

		// Number of comb filters in parallel network
		var numCombs = 4;

		// Master modulator
		var mod = LFTri.kr(modfreq);

		// Feedback in
		fb = LocalIn.ar(numCombs) * fbamount.linlin(0.0,1.0,0.0,0.25);

		// Parallel comb filter network (schroeder inspired)
		sig = Array.fill(numCombs, {|num|

			var thisfb = LPF.ar(fb[num], damp.linexp(0.0,1.0,20000.0,{rrand(0.9,1.1)} * 200.0));

			CombC.ar(
				thisfb + source, 
				maxdelaytime: 0.2, 
				delaytime: SinOsc.kr(
					num+1/(100) * {rrand(0.90000234, 1.1052)} * mod.linlin(-1.0,1.0,0.0,8.0), 
					2pi.rrand(-2pi)
				).range(1 - modamount,1.0 + modamount) // modulation amount here
				* {rrand(0.01,0.09)}
				* verbtime,
				decaytime: 1.0 + verbtime
			)

		}).sum / numCombs;

		// Add diffusion
		sig = AllpassC.ar(sig, 0.2, 0.017123, 1.0);
		sig = AllpassC.ar(sig, 0.2, 0.0512385, 1.0);

		// Feedback out
		LocalOut.ar(
			LeakDC.ar(
				AllpassC.ar(
					sig, 
					0.2, 
					(SinOsc.kr((modfreq+fbamount) / 100 * LFNoise2.kr(1).range(0.5,1.0) * mod * verbtime.reciprocal).range(0.01,0.025) * modamount).clip(0.0001,10.0),
					decaytime: verbtime + 0.5
				).tanh
			)
		);

		// Filter output
		sig = LPF.ar(sig, highcut);
		sig = HPF.ar(sig, lowcut);

		^sig
	}
}


// Stereo reverb
// The starting point for this one can be found here: https://en.wikibooks.org/wiki/Designing_Sound_in_SuperCollider/Schroeder_reverb
Skyrverb2{
	*ar{|in, verbtime=2.31, modfreq=0.05, damp=0.0008, decay=0.35, lowcut=80, highcut=8500, modamount=0.25|	
		var input, output, delrd, sig, deltimes, modulators;

		verbtime = verbtime.lag;

		input = in;

		// Read our 4-channel delayed signals back from the feedback loop
		delrd = LPF.ar(LeakDC.ar(LocalIn.ar(4) * decay).tanh, damp.linexp(0.0,1.0,20000.0,50.0));

		// This will be our eventual output, which will also be recirculated
		output = input + delrd[[0,1]];
		output = LPF.ar(output, damp.linexp(0.0,1.0,20000.0,50.0));

		output = AllpassC.ar(output, 0.2, 0.015808725182 * verbtime );
		output = AllpassC.ar(output, 0.2, 0.035808725182 * verbtime );

		// Cross-fertilise the four delay lines with each other:
		sig = [output[0]+output[1], output[0]-output[1], delrd[2]+delrd[3], delrd[2]-delrd[3]];
		sig = [sig[0]+sig[2], sig[1]+sig[3], sig[0]-sig[2], sig[1]-sig[3]];

		modulators = SinOsc.ar([modfreq, modfreq*1.71719, modfreq*0.98571234, modfreq*0.81818], {2pi.rrand(-2pi)}!4);

		// Attenutate the delayed signals so they decay:
		sig = sig * [0.4, 0.37, 0.333, 0.3];

		// Here we give delay times in milliseconds, convert to seconds,
		// then compensate with ControlDur for the one-block delay
		// which is always introduced when using the LocalIn/Out fdbk loop
		deltimes = [101, 143, 165, 177] * 0.001 - ControlDur.ir * verbtime;

		// Modulation
		deltimes = deltimes * modulators.range(1-modamount,1+modamount); 

		// Apply the delays and send the signals into the feedback loop
		LocalOut.ar(DelayC.ar(sig, deltimes, deltimes));

		// Now let's hear it:
		output = output / 4;
		output = LPF.ar(output, highcut);
		output = HPF.ar(output, lowcut);

		^output
	}
}
