# Polyhedron Generation Tool for unity

![Screenshot](ScreenShot.png)

Provides function to generate simple, and procedurally modify existing polhedron shapes.
Can store these shapes as both Unity Meshobjects and as custom face and neighbor based scriptable Object.
Operations include geometic trucation, dual, spherize, and for triangluar faces- tesselation.

## Dependencies

  Requires Unity
  Packages:
	UniTask by Cysharp https://github.com/Cysharp/UniTask.git?path=src/UniTask/Assets/Plugins/UniTask
    	CommonAssetTypes https://github.com/glurth/CommonAssetTypes.git

## Installation

In the unity editor, package manager window, click to add a package from a git url. 
First add the dependencies specified in this file, then add:https://github.com/glurth/UIPrefabGenerator.git

## Usage

One may use the polyhedron and associated sample assets directly, or may use the interace in the main scene, to generate custom ones.

The scene provides an interface that will allow the user to perform operations on a "starting polyhedron".  Currently, the user will need to specify the starting object, before running the scene, in the "Maker" scene object's inspector.  User may specify a starting Polyhedron asset, or a Mesh asset that the system will attempt to generate a polyhedron from (or neither to use the default starting point of an icosahedron).


## File structure

RunTime:
	Contains the code files and class definitions.  
SampleAssets:
 	Contains a bunch of pregenerated unity and polyhedron assets.  Each Polyhedron contains the data that defines a polyhedron.  It MAY also store this information in two other forms, for efficiency: a regular unity mesh, and a FacesAndNeighbors asset which proves usful for path finding and maps.
	Two materials with thier own custom shaders, and a test texture.  These materials can been tuned to show how both UV0 and UV1 are mapped.
	The main scene with an interface to create new polyhedrons.
	The polyhedrons contain in here started with an icosahedron and used multiple triangular face teselation operations, and a final dual operation when generated.
	The TruncPyramidFolder contains polyhedrons generated, starting with an icosahedron, using the truncate corner, and pyramid face operations in an alternating sequence.

## License

  This code is proprietary, and no license is granted without written permission. Free Licence for my fellow indie devs, just email me at glurth at gmail dot com.


#Contributing

  Contributions are welcome! Please submit pull requests or open issues for any bugs or enhancements.
Contact

For any questions or feedback, feel free to reach out.