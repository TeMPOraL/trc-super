-module(ray_scene).
-compile(export_all).

%%%% Raycaster Scene Objects.
%%%% This implementation supports following objects:
%%%% Object category:
%%%%	-> Sphere
%%%%	-> Plane
%%%% Light source category:
%%%%	-> Point light source

%%% Point-in-object tests
%% Test if point is inside a sphere object.
point_inside([sphere, _, InvTransform, Radius | _], Vector) ->
	{Tx, Ty, Tz} = ray_math:matrix4_transform_vector3_point(InvTransform, Vector),
	case Tx*Tx + Ty*Ty + Tz*Tz of
		Val when Val < Radius*Radius ->
			true;
		Val -> false
	end;
%% Test if point is inside a plane object - of course it never happens.
point_inside([plane | _], Vector) -> false.

%%% Ray-Object intersections
%% Get intersections for ray and sphere.
intersections([sphere, Transform, InvTransform, Radius | _], {RayStart, RayDir}) ->
	RayStartTr = ray_math:matrix4_transform_vector3_point(InvTransform, RayStart),
	RayDirTr = ray_math:matrix4_transform_vector3_direction(InvTransform, RayDir),
	%%....