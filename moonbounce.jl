### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ e910b348-9e5a-40d9-be9a-e78aacc47743
using Pkg; Pkg.activate("."); Pkg.instantiate()

# ╔═╡ feb48b54-7423-11ec-2426-23da81f3c079
using SatelliteToolbox, Dates

# ╔═╡ c3e74733-4f49-4849-8d64-5def3f5392e9
using LinearAlgebra,Roots

# ╔═╡ 1bd4c734-e245-4027-ac70-2d8043922cbe
using Plots, Plots.PlotMeasures

# ╔═╡ 274cb358-741c-4b43-bd6a-f96f3453c7fa
md"Download some Earth orientation parameters:"

# ╔═╡ 5c813a37-6916-4214-bb0a-380ca33a767f
eop_IAU1980 = get_iers_eop();

# ╔═╡ 5568e66d-8d62-4b78-825a-ec9e382752ea
md"### Some Helpful Definitions"

# ╔═╡ a8e6acad-212e-4cbe-b2b9-03b7292636f0
"""
Speed of light in a vacuum
"""
c = 299792458 # m/s

# ╔═╡ cedb1362-47a2-44ed-abd8-34b37c9d009e
"""
The frequency of the LoRa signal used in the experiment was between 430 and 440 MHz
"""
fLoRa = 435e6 # MHz

# ╔═╡ e271d16d-d600-41c5-9d58-3b56d1683c04
"""
Approximate latitude, longitue, and altitude of the Dwingeloo Radio Telescope used to transmit the signal
"""
dwingelooLLA = (52.8121, 6.3971, 0) # 0 meters elevation

# ╔═╡ 9b2af5c2-42ea-4b4f-a2e7-9b67a98da790
"""
    earthLLA_r(latdeg, lngdeg, altkm, jd)

J2000 position on Julian day `jd` of a point on the Earth's surface specified by a latitude, longitude, and altitude above the geoid in kilometers.
"""
function earthLLA_r(latdeg, lngdeg, altkm, jd)

	rITRFtoJ2K = r_ecef_to_eci(ITRF(), J2000(), jd, eop_IAU1980)
	itrfPos = geodetic_to_ecef(deg2rad(latdeg), deg2rad(lngdeg), altkm*1000)

	posJ2K = rITRFtoJ2K * itrfPos
end

# ╔═╡ b4fc572f-a9e8-4dac-b7b8-023a4229429d
"""
    earthLLA_rv(latdeg, lngdeg, altkm, jd)

Position and velocity of a point on the Earth's surface on the Julian day `jd` in the J2000 frame. Velocity is derived from a forward finite difference of position.
"""
function earthLLA_rv(latdeg, lngdeg, altkm, jd)
	dt = 1 # seconds
	pos1 = earthLLA_r(latdeg, lngdeg, altkm, jd)
	pos2 = earthLLA_r(latdeg, lngdeg, altkm, jd+dt/86400)

	pos1, (pos2-pos1)/dt
end

# ╔═╡ 11e63381-a19c-4202-96bb-3588c9d8d3a5
"""
    moon_rv(jd)

Position and velocity of the Moon at the Julian date `jd` in the J2000 frame (ignoring conversion between UTC time and Barycentric Dynamical Time). Velocity
is derived from a forward finite difference of position.
"""
function moon_rv(jd)
	R = r_eci_to_eci(MOD(), J2000(), jd)

	dt = 1 # seconds
	
	jd2 = jd + dt/86400

	# convert MOD positions to J2000
	# TODO: convert JD_UTC to JD_TBD
	pos =  R*moon_position_i(jd)
	pos2 = R*moon_position_i(jd2)

	pos, (pos2-pos)/dt
end

# ╔═╡ 57fc4f13-2218-43fe-8193-350f6b6db12e
"""

    moonLLA_r(latdeg, lngdeg, altkm, jd)

J2000 position on Julian day `jd` of a point on the Moon's surface specified by a latitude, longitude, and an altitude above the Moon's average radius in kilometers.

I've defined the Moon frame has +x pointing from Moon center to Earth, +z pointing in the Moon's orbital angular momentum direction, and +y completing the right-handed set.
"""
function moonLLA_r(latdeg, lngdeg, altkm, jd)
	rMoon = 1738.1e3 # meters

	# assumes moon rv is in J2000
	moonpos, moonvel = moon_rv(jd)

	xhat = -normalize(moonpos)
	zhat = normalize(cross(moonpos, moonvel))
	yhat = normalize(cross(zhat, xhat))

	rMoonToJ2000 = [xhat yhat zhat]

	# convert polar angles to a Cartesian position in the Moon frame
	posMoonFrame = (rMoon + altkm*1000) * [cosd(latdeg)*cosd(lngdeg),
										   cosd(latdeg)*sind(lngdeg),
										   sind(latdeg)]
	# conver from my Moon frame to J2000
	posJ2000 = rMoonToJ2000 * posMoonFrame + moonpos
end

# ╔═╡ 73702699-2464-4b22-abc9-422a2285bd4d
"""
    moonLLA_rv(latdeg, lngdeg, altkm, jd)

Position and velocity of a point on the Moons's surface on the Julian day `jd` in the J2000 frame ignoring Barycentric Dynamical Time conversion. Velocity is derived from a forward finite difference of position.
"""
function moonLLA_rv(latdeg, lngdeg, altkm, jd)
	dt = 1 # second
	pos1 = moonLLA_r(latdeg, lngdeg, altkm, jd)
	pos2 = moonLLA_r(latdeg, lngdeg, altkm, jd+dt/86400)

	pos1, (pos2-pos1)/dt
end

# ╔═╡ 461c3f7b-c2a0-4ae2-849b-026978f4dfe3
md"""
## The Main Event: `moonbounce(...)`

"""

# ╔═╡ 0e00242a-1023-480f-a187-c5e601712d59
function moonbounce(earthlat, earthlng, earthaltkm,  moonlat, moonlng, moonaltkm,  transmitjd; assumedLightDelaySeconds = 1.0)

	transmitPos, transmitVel = earthLLA_rv(earthlat, earthlng, earthaltkm, transmitjd)
	
	function earthToMoonDelay(delaySeconds)
		moonPos, moonVel = moonLLA_rv(moonlat, moonlng, moonaltkm, transmitjd + delaySeconds/86400)

		dist = norm(moonPos - transmitPos)

		dist/c - delaySeconds
	end

	e2mSeconds = find_zero(earthToMoonDelay, assumedLightDelaySeconds)

	# time of signal reflection by a reflector on the Moon
	# and position and velocity of the reflector at that time
	reflectionTime = transmitjd + e2mSeconds/86400
	reflectPos, reflectVel = moonLLA_rv(moonlat, moonlng, moonaltkm, reflectionTime)

	function moonToEarthDelay(delaySeconds)
		earthPos, earthVel = earthLLA_rv(earthlat, earthlng, earthaltkm, reflectionTime + delaySeconds/86400)

		dist = norm(earthPos - reflectPos)

		dist/c - delaySeconds
	end

	m2eSeconds = find_zero(moonToEarthDelay, assumedLightDelaySeconds)

	# time of reception back on Earth,
	# and the position and velocity of the receiver at that time
	receptionTime = reflectionTime + m2eSeconds/86400
	receptionPos, receptionVel = earthLLA_rv(earthlat, earthlng, earthaltkm, receptionTime)

	(transmission=(time=transmitjd, pos=transmitPos, vel=transmitVel),
	 reflection=(time=reflectionTime, pos=reflectPos,  vel=reflectVel),
	 reception= (time=receptionTime,  pos=receptionPos,vel=receptionVel))
end

# ╔═╡ f3651d01-81bf-4e3e-ab60-cd70541ce24d
"""
    function dopplerFactorBetween(a, b)
Calculates the Doppler effect frequency shift factor of a light wave transmitted from one point and received by another, both with instantaneous velocities.

Based on the equation from the Wikipedia page [Doppler effect](https://en.wikipedia.org/wiki/Doppler_effect#General).
"""
function dopplerFactorBetween(a, b)
	# vector between the two points
	rAtoB = b.pos-a.pos
	rHat = normalize(rAtoB)

	# source velocity relative to the receiver
	vS = dot(a.vel, rHat)

	# receiver velocity, relative to source
	vR = dot(b.vel, -rHat)  # -rHat because we want the direction from b to a

	factor = (c+vR)/(c-vS)
end

# ╔═╡ 5b4fd441-452e-4bf2-b606-b4724909cdd7
md"
## Generating Points on the Moon

I'd like to evenly cover the Moon with points that will act as reflector locations. A great way to evenly distribute a number of points around a sphere is to use, of all things, the Fibonacci spiral. I'm following the procedure described by [extremelearning.com.au](http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/)

"

# ╔═╡ e019b7ef-3671-4c48-97a0-aed0c0857efc
n = 500

# ╔═╡ c6af0b6f-955f-4ab3-a105-4996aa1a5383
# generate some reflector latitudes and longitudes
refLats, refLngs = let
	ϕ = (1 + 5^0.5)/2
	i = 0:n-1

	lng = rad2deg.(2pi/ϕ * i)
	# colatitude - zero is +z (the north pole), pi is -z (south pole)
	colat = acosd.(1 .- 2/n .* (i.+0.5))

	colat, lng
end

# ╔═╡ 303254ab-7ab0-476e-892c-a2866a48db29
md"""
## Simulate a bounce

With the `moonbounce` function in place, all that needs to be done is turn each of the transmit-reflect-receive sets it produces into a delay-Doppler pair. I'm also returning a "strength" result for each bounce, which hopefully is an approximate indication of how strong a reflection each bounce is. Bouncing off a part of the Moon facing directly at the transmitter should give a strong reflection, while a glancing bounce should be harder to see in the final plot.

Also included here is the handy `angleBetweend` function that returns the angle between two 3D vectors in degrees.
"""

# ╔═╡ 80c79503-c342-4f31-9356-bae9a29aaf4f
"""
    angleBetweend(a,b)
Angle between two (unnormalized) vectors in degrees

Examples:
```jldoctest
julia> using LinearAlgebra

julia> angleBetweend([1,0,0], [0,1,0])
90.0

julia> angleBetweend([1,0,0], [1,1,0])
45.0
```
"""
function angleBetweend(a,b)
	na = normalize(a)
	nb = normalize(b)
	atand(norm(cross(na, nb)), dot(na,nb))
end;

# ╔═╡ e165d46b-af7f-47aa-8fcb-5999606bdb0f
function dopplerDelayPairsAt(jd)

	# generate a Doppler-delay pair a signal transmitted at time t bouncing of each
	# reflector on the Moon
	dopplerDelayPairs = map(refLats, refLngs) do moonLat, moonLng
		bounce = moonbounce(dwingelooLLA..., moonLat,moonLng,0, jd)
	
		moonPosAtReflection = moon_rv(bounce.reflection.time)[1]
		reflectionEarthElevation = 90-angleBetweend(bounce.reflection.pos-moonPosAtReflection, bounce.transmission.pos-bounce.reflection.pos)
	
		transmitterMoonElevation = 90-angleBetweend(bounce.reflection.pos-bounce.transmission.pos, bounce.transmission.pos)
	
		if(reflectionEarthElevation > 5 && # skip "reflections" that have a low (or
										   # negative!) elevation on the Moon
			transmitterMoonElevation>10)   # and also skip ones that have too low 
										   # an elevation on the Earth
			
			shift = fLoRa * 
					(dopplerFactorBetween(bounce.transmission, bounce.reflection) *
					 dopplerFactorBetween(bounce.reflection, bounce.reception)
					 - 1)

			# convert delay in Julian days to delay in seconds
			delay = (bounce.reception.time - bounce.transmission.time)*86400
	
			shift, delay, sind(reflectionEarthElevation)
		else
			nothing, nothing, nothing
		end
	end

	# filter out all the passes that had an elevation that was too low
	dopplers = filter(x->!isnothing(x), map(x->x[1], dopplerDelayPairs));
	delays = filter(x->!isnothing(x), map(x->x[2], dopplerDelayPairs));
	strengths = filter(x->!isnothing(x), map(x->x[3], dopplerDelayPairs));

	dopplers, delays, strengths
end

# ╔═╡ 21e56be6-dafc-46e6-9831-b258d2d7aedb
function dopplerDelayPlotAt(jd; plotargs...)
	dopplers, delays, strengths = dopplerDelayPairsAt(jd)
	
	title = Dates.format(DateTime(jd_to_date(jd)[1:end-1]...), "u. dd HH:MM")
	
	scatter(dopplers, delays, 
		alpha=strengths,     # first cut at Lambertian-type reflection strength TODO
		label=false,
		xlabel="Frequency shift (Hz)",
		ylabel="Delay since start of transmission (s)",
		markerstrokewidth=0,
		guidefont = 8, tickfont = 8,
		title = title;
		plotargs...
	)
end

# ╔═╡ ece89996-7e66-4923-81d2-549bd332752e
let
	transmitTimes = (DateTime(2021, 10, 5, 3, 24, 55):Dates.Hour(1):DateTime(2021, 10, 5, 18, 24, 55))
	
	bouncePlots = dopplerDelayPlotAt.(date_to_jd.(transmitTimes), aspect_ratio=400, ylabel="Delay (s)", xrot=40)
	
	p = plot(bouncePlots..., layout=(4,4),size=(1400,1000), bottom_margin=20px, left_margin=20px)

	titleplot = plot(title="Moon Bounce Doppler Delay October 5, 2021", grid=false, showaxis=false, xaxis=nothing, yaxis=nothing, bottom_margin=-20Plots.px)
	
	plot(titleplot, p, layout=@layout([A{0.01h}; B]))
end

# ╔═╡ 40f45a40-01d2-4569-8825-63a81779cd99
md"
## A Plot that Matches the Original

By guessing and checking, I found that a transimission time of 12:24:55 UTC on October 5, 2021 produces a plot that is a pretty good match to the one shared by Lacuna [here](https://lacuna.space/lora-moon-bounce/)
"


# ╔═╡ 59550d8b-13de-4d99-9164-b5b0a42fea7e
t = date_to_jd(DateTime(2021, 10, 5, 12, 24, 55))

# ╔═╡ 833c3151-9e37-4f72-9ffa-bf1958d9adc7
bounce = moonbounce(dwingelooLLA..., 0,0,0, t)

# ╔═╡ 7059f0f2-c620-4353-a2c7-2070d90320fe
dopplerFactorBetween(bounce.transmission, bounce.reflection)

# ╔═╡ 494659ca-1258-4f04-a8b3-c7889610ab43
dopplerFactorBetween(bounce.reflection, bounce.reception)

# ╔═╡ 356a3352-4a77-41e2-9206-749ca8d40983
# take two Doppler factors and convert to a frequncy shift in Hertz
shift = fLoRa * (dopplerFactorBetween(bounce.transmission, bounce.reflection) *dopplerFactorBetween(bounce.reflection, bounce.reception) - 1)

# ╔═╡ ce1238ea-0897-4d2b-902a-14da86387833
delay = bounce.reception.time - bounce.transmission.time

# ╔═╡ fd28654f-18aa-4e27-a83d-91f2e6cd81c4
dopplerDelayPlotAt(t)

# ╔═╡ Cell order:
# ╠═e910b348-9e5a-40d9-be9a-e78aacc47743
# ╠═feb48b54-7423-11ec-2426-23da81f3c079
# ╠═c3e74733-4f49-4849-8d64-5def3f5392e9
# ╠═1bd4c734-e245-4027-ac70-2d8043922cbe
# ╟─274cb358-741c-4b43-bd6a-f96f3453c7fa
# ╠═5c813a37-6916-4214-bb0a-380ca33a767f
# ╟─5568e66d-8d62-4b78-825a-ec9e382752ea
# ╠═a8e6acad-212e-4cbe-b2b9-03b7292636f0
# ╠═cedb1362-47a2-44ed-abd8-34b37c9d009e
# ╠═e271d16d-d600-41c5-9d58-3b56d1683c04
# ╠═9b2af5c2-42ea-4b4f-a2e7-9b67a98da790
# ╠═b4fc572f-a9e8-4dac-b7b8-023a4229429d
# ╠═11e63381-a19c-4202-96bb-3588c9d8d3a5
# ╠═57fc4f13-2218-43fe-8193-350f6b6db12e
# ╠═73702699-2464-4b22-abc9-422a2285bd4d
# ╠═461c3f7b-c2a0-4ae2-849b-026978f4dfe3
# ╠═0e00242a-1023-480f-a187-c5e601712d59
# ╠═f3651d01-81bf-4e3e-ab60-cd70541ce24d
# ╠═833c3151-9e37-4f72-9ffa-bf1958d9adc7
# ╠═7059f0f2-c620-4353-a2c7-2070d90320fe
# ╠═494659ca-1258-4f04-a8b3-c7889610ab43
# ╠═356a3352-4a77-41e2-9206-749ca8d40983
# ╠═ce1238ea-0897-4d2b-902a-14da86387833
# ╟─5b4fd441-452e-4bf2-b606-b4724909cdd7
# ╠═e019b7ef-3671-4c48-97a0-aed0c0857efc
# ╠═c6af0b6f-955f-4ab3-a105-4996aa1a5383
# ╟─303254ab-7ab0-476e-892c-a2866a48db29
# ╠═80c79503-c342-4f31-9356-bae9a29aaf4f
# ╠═e165d46b-af7f-47aa-8fcb-5999606bdb0f
# ╠═21e56be6-dafc-46e6-9831-b258d2d7aedb
# ╠═ece89996-7e66-4923-81d2-549bd332752e
# ╟─40f45a40-01d2-4569-8825-63a81779cd99
# ╠═59550d8b-13de-4d99-9164-b5b0a42fea7e
# ╠═fd28654f-18aa-4e27-a83d-91f2e6cd81c4
