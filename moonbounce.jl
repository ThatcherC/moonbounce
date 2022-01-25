### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ feb48b54-7423-11ec-2426-23da81f3c079
using SatelliteToolbox, Dates

# ╔═╡ c3e74733-4f49-4849-8d64-5def3f5392e9
using LinearAlgebra

# ╔═╡ 7e991914-1050-4e85-b729-284f7f68e457
using Roots

# ╔═╡ 9005eda7-eed6-487e-9033-25780cfccb79


# ╔═╡ 5c813a37-6916-4214-bb0a-380ca33a767f
eop_IAU1980 = get_iers_eop();

# ╔═╡ a8e6acad-212e-4cbe-b2b9-03b7292636f0
# lightspeed
c = 299792458 # m/s

# ╔═╡ cedb1362-47a2-44ed-abd8-34b37c9d009e
fLoRa = 435e6 # Hz, LoRa operated at 430-440 MHz in the LoRa moonbounce experiment

# ╔═╡ 11e63381-a19c-4202-96bb-3588c9d8d3a5
function moon_rv(date::DateTime)
	# TODO: convert Moon coords to J2000
	R = I

	jd = date_to_jd(date)
	jd2 = date_to_jd(date+Dates.Second(1))
	
	pos = R*moon_position_i(jd)
	pos2 = R*moon_position_i(jd2)

	pos, pos2-pos # assumes dt is 1 second
end

# ╔═╡ 9b2af5c2-42ea-4b4f-a2e7-9b67a98da790
function earthLLA_r(latdeg, lngdeg, altkm, date::DateTime)

	rITRFtoJ2K = r_ecef_to_eci(ITRF(), J2000(), date_to_jd(date), eop_IAU1980)
	itrfPos = geodetic_to_ecef(deg2rad(latdeg), deg2rad(lngdeg), altkm*1000)

	posJ2K = rITRFtoJ2K * itrfPos
end

# ╔═╡ b4fc572f-a9e8-4dac-b7b8-023a4229429d
function earthLLA_rv(latdeg, lngdeg, altkm, date::DateTime)
	pos1 = earthLLA_r(latdeg, lngdeg, altkm, date)
	pos2 = earthLLA_r(latdeg, lngdeg, altkm, date+Dates.Second(1))

	pos1, pos2-pos1 # assumes dt is 1 second
end

# ╔═╡ 57fc4f13-2218-43fe-8193-350f6b6db12e
"""

Moon frame has x pointing from moon center to Earth, y pointing in negative moon velocity direction
"""
function moonLLA_r(latdeg, lngdeg, altkm, date::DateTime)
	rMoon = 1738.1e3 # meters

	# assumes moon rv is in J2000
	moonpos, moonvel = moon_rv(date)

	xhat = -normalize(moonpos)
	zhat = normalize(cross(moonpos, moonvel))
	yhat = normalize(cross(zhat, xhat))

	rMoonToJ2000 = [xhat yhat zhat]

	posMoonFrame = (rMoon + altkm*1000) * [cosd(latdeg)*cosd(lngdeg),
										   cosd(latdeg)*sind(lngdeg),
										   sind(latdeg)]
	posJ2000 = rMoonToJ2000 * posMoonFrame + moonpos
end

# ╔═╡ 73702699-2464-4b22-abc9-422a2285bd4d
function moonLLA_rv(latdeg, lngdeg, altkm, date::DateTime)
	pos1 = moonLLA_r(latdeg, lngdeg, altkm, date)
	pos2 = moonLLA_r(latdeg, lngdeg, altkm, date+Dates.Second(1))

	pos1, pos2-pos1
end

# ╔═╡ 0e00242a-1023-480f-a187-c5e601712d59
function moonbounce(earthlat, earthlng, earthaltkm,  moonlat, moonlng, moonaltkm,  transmittime::DateTime; assumedLightDelayMs = 1000)

	transmitPos, transmitVel = earthLLA_rv(earthlat, earthlng, earthaltkm, transmittime)
	
	function earthToMoonDelayMillis(delaymillis)
		moonPos, moonVel = moonLLA_rv(moonlat, moonlng, moonaltkm, transmittime + Dates.Millisecond(round(delaymillis)))

		dist = norm(moonPos - transmitPos)

		dist/c*1000 - delaymillis
	end

	e2mMillis = find_zero(earthToMoonDelayMillis, assumedLightDelayMs)

	# time of signal reflection by a reflector on the Moon
	# and position and velocity of the reflector at that time
	reflectionTime = transmittime + Dates.Millisecond(round(e2mMillis))
	reflectPos, reflectVel = moonLLA_rv(moonlat, moonlng, moonaltkm, reflectionTime)

	function moonToEarthDelayMillis(delaymillis)
		earthPos, earthVel = earthLLA_rv(earthlat, earthlng, earthaltkm, reflectionTime + Dates.Millisecond(round(delaymillis)))

		dist = norm(earthPos - reflectPos)

		dist/c*1000 - delaymillis
	end

	m2eMillis = find_zero(moonToEarthDelayMillis, assumedLightDelayMs)

	# time of reception back on Earth,
	# and the position and velocity of the receiver at that time
	receptionTime = reflectionTime + Dates.Millisecond(round(m2eMillis))
	receptionPos, receptionVel = earthLLA_rv(earthlat, earthlng, earthaltkm, receptionTime)

	(transmission=(time=transmittime, pos=transmitPos, vel=transmitVel),
	 reflection=(time=reflectionTime, pos=reflectPos,  vel=reflectVel),
	 reception= (time=receptionTime,  pos=receptionPos,vel=receptionVel))
end

# ╔═╡ 59550d8b-13de-4d99-9164-b5b0a42fea7e
t = DateTime(2021, 10, 24, 10, 24)

# ╔═╡ e271d16d-d600-41c5-9d58-3b56d1683c04
# approximate location of the Dwingeloo Radio Observatory
dwingelooLLA = (52.8121050606517, 6.3971174915299684, 0) # 0 meters elevation

# ╔═╡ f3651d01-81bf-4e3e-ab60-cd70541ce24d
"""
Calculates the Doppler effect frequency shift factor of a light wave transmitted from one point and received by another, both with instantaneous velocities.

https://en.wikipedia.org/wiki/Doppler_effect#General
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

# ╔═╡ 833c3151-9e37-4f72-9ffa-bf1958d9adc7
bounce = moonbounce(dwingelooLLA..., 0,0,0, t+Dates.Hour(-8))

# ╔═╡ 7059f0f2-c620-4353-a2c7-2070d90320fe
dopplerFactorBetween(bounce.transmission, bounce.reflection)

# ╔═╡ 494659ca-1258-4f04-a8b3-c7889610ab43
dopplerFactorBetween(bounce.reflection, bounce.reception)

# ╔═╡ 356a3352-4a77-41e2-9206-749ca8d40983
shift = fLoRa * (dopplerFactorBetween(bounce.transmission, bounce.reflection) *dopplerFactorBetween(bounce.reflection, bounce.reception) - 1)

# ╔═╡ ce1238ea-0897-4d2b-902a-14da86387833
delay = bounce.reception.time - bounce.transmission.time

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SatelliteToolbox = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"

[compat]
Roots = "~1.3.14"
SatelliteToolbox = "~0.9.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "54fc4400de6e5c3e27be6047da2ef6ba355511f8"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.6"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "67551df041955cc6ee2ed098718c8fcd7fc7aebe"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.12.0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OptionalData]]
git-tree-sha1 = "d047cc114023e12292533bb822b45c23cb51d310"
uuid = "fbd9d27c-2d1c-5c1c-99f2-7497d746985d"
version = "1.0.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.ReferenceFrameRotations]]
deps = ["Crayons", "LinearAlgebra", "Printf", "StaticArrays"]
git-tree-sha1 = "d526371cec370888f485756a4bf8284ab531860b"
uuid = "74f56ac7-18b3-5285-802d-d4bd4f104033"
version = "1.0.1"

[[deps.RemoteFiles]]
deps = ["Dates", "FileIO", "HTTP"]
git-tree-sha1 = "54527375d877a64c55190fb762d584f927d6d7c3"
uuid = "cbe49d4c-5af1-5b60-bb70-0a60aa018e1b"
version = "0.4.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "0abe7fc220977da88ad86d339335a4517944fea2"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.14"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SatelliteToolbox]]
deps = ["Crayons", "Dates", "DelimitedFiles", "Interpolations", "LinearAlgebra", "OptionalData", "Parameters", "PolynomialRoots", "PrettyTables", "Printf", "Reexport", "ReferenceFrameRotations", "RemoteFiles", "SparseArrays", "StaticArrays", "Statistics"]
git-tree-sha1 = "1831cced8785398bf38577e8cf46380d349cf4c9"
uuid = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
version = "0.9.4"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2884859916598f974858ff01df7dfc6c708dd895"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═9005eda7-eed6-487e-9033-25780cfccb79
# ╠═feb48b54-7423-11ec-2426-23da81f3c079
# ╠═c3e74733-4f49-4849-8d64-5def3f5392e9
# ╠═7e991914-1050-4e85-b729-284f7f68e457
# ╠═5c813a37-6916-4214-bb0a-380ca33a767f
# ╠═a8e6acad-212e-4cbe-b2b9-03b7292636f0
# ╠═cedb1362-47a2-44ed-abd8-34b37c9d009e
# ╠═11e63381-a19c-4202-96bb-3588c9d8d3a5
# ╠═9b2af5c2-42ea-4b4f-a2e7-9b67a98da790
# ╠═b4fc572f-a9e8-4dac-b7b8-023a4229429d
# ╠═57fc4f13-2218-43fe-8193-350f6b6db12e
# ╠═73702699-2464-4b22-abc9-422a2285bd4d
# ╠═0e00242a-1023-480f-a187-c5e601712d59
# ╠═59550d8b-13de-4d99-9164-b5b0a42fea7e
# ╠═e271d16d-d600-41c5-9d58-3b56d1683c04
# ╠═f3651d01-81bf-4e3e-ab60-cd70541ce24d
# ╠═833c3151-9e37-4f72-9ffa-bf1958d9adc7
# ╠═7059f0f2-c620-4353-a2c7-2070d90320fe
# ╠═494659ca-1258-4f04-a8b3-c7889610ab43
# ╠═356a3352-4a77-41e2-9206-749ca8d40983
# ╠═ce1238ea-0897-4d2b-902a-14da86387833
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
