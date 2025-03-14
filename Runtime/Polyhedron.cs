using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Cysharp.Threading.Tasks;


class GeometryException : System.Exception
{
    public GeometryException(string message) : base(message)
    {
    }
}

/// <summary>
/// Stores mesh information in a class that can be instantiated and used outside of the main unity thread.
/// </summary>
public class MeshData
{
    public UnityEngine.Rendering.IndexFormat indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
    public Vector3[] meshVerts;
    public int[] meshTris;
    public Vector3[] meshNormals;

    public Vector2[] meshUV0s;
    public Vector2[] meshUV1s;
    public Vector2[] meshUV2s;
    public FacesAndNeighbors facesAndNeighborsRef;
    public string name;
    public Mesh ToMesh()
    {
        /*
            SphereizeMeshUV0s(meshToUseRef);
        meshToUseRef.name = faces.Count.ToString() + "FacePoly";
            if (doFacesAndNeighbors)
                facesAndNeighbors.meshRef = meshToUseRef;
        */
        Mesh newMesh = new Mesh();
        newMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        newMesh.SetVertices(meshVerts);
        newMesh.SetTriangles(meshTris, 0);
        newMesh.SetNormals(meshNormals);
        newMesh.SetUVs(0, meshUV0s);
        newMesh.SetUVs(1, meshUV1s);
        newMesh.SetUVs(2, meshUV2s);
        
        newMesh.name = name;
        if (facesAndNeighborsRef!=null)
            facesAndNeighborsRef.meshRef = newMesh;

        return newMesh;
    }
}


public static class Vector3Extensions
{
    public static Vector3 AvgPos(this List<Vector3> points)
    {
        int count =0;
        Vector3 vectorSum= Vector3.zero;
        foreach (Vector3 p in points)
        {
            vectorSum += p;
            count++;
        }
        return vectorSum / count;
    }

    /// <summary>
    /// Determines whether two floating-point numbers are approximately equal based on a fractional difference.
    /// This extension method is useful for comparing floating-point numbers where precision errors may occur.
    /// </summary>
    /// <param name="a">The first floating-point number.</param>
    /// <param name="b">The second floating-point number.</param>
    /// <param name="fractionalDiff">The maximum fractional difference allowed for the numbers to be considered equal. Default is 0.001.</param>
    /// <returns>True if the numbers are approximately equal; otherwise, false.</returns>
    public static bool CloseEqual(this float a, float b, float fractionalDiff = .001f)
    {
        float aMinusB = a - b;
        if (aMinusB < 0)
            aMinusB *= -1; // Convert the difference to absolute value
        if (a < 0)
            a *= -1; // Convert 'a' to its absolute value
        if (b < 0)
            b *= -1; // Convert 'b' to its absolute value

        float larger = a;
        if (a < b) larger = b; // Determine the larger of the two values

        // Use the larger value as a reference for the fractional difference to ensure the comparison is scale-invariant.
        if (aMinusB <= fractionalDiff * larger)
            return true;
        return false;
    }

    /// <summary>
    /// Determines whether two Vector3 values are approximately equal based on a fractional difference.
    /// This extension method compares each component (x, y, z) of the vectors using the CloseEqual method.
    /// </summary>
    /// <param name="a">The first Vector3 value.</param>
    /// <param name="b">The second Vector3 value.</param>
    /// <param name="fractionalDiff">The maximum fractional difference allowed for each component to be considered equal. Default is 0.001.</param>
    /// <returns>True if all components are approximately equal; otherwise, false.</returns>
    public static bool CloseEqual(this Vector3 a, Vector3 b, float fractionalDiff = .001f)
    {
        if (!CloseEqual(a.x, b.x, fractionalDiff)) return false;
        if (!CloseEqual(a.y, b.y, fractionalDiff)) return false;
        if (!CloseEqual(a.z, b.z, fractionalDiff)) return false;
        return true;
    }


    /// <summary>
    /// Returns the circular index of the given integer `a` within the specified range of 0 to `size - 1`.
    /// If `a` exceeds or goes below the bounds, it wraps around to ensure the result is always between 0 and `size - 1`.
    /// For example:
    /// - If `a` is greater than or equal to `size`, it will wrap around to the beginning (0).
    /// - If `a` is negative, it will wrap around to the end (`size - 1`).
    /// </summary>
    /// <param name="a">The index to be wrapped around.</param>
    /// <param name="size">The size of the range (exclusive upper bound).</param>
    /// <returns>The circularly wrapped index.</returns>
    public static int CircularIndex(this int a, int size)
    {
        int b = a % size;
        if (b < 0) b += size;
        return b;
    }

    const float radToRots = 1.0f / (2f * Mathf.PI);
    /// <summary>
    /// Converts a 3D vertex on the unit sphere to its corresponding cylindrical UV coordinates.
    /// The longitude corresponds to the rotation around the Y-axis (north pole), and the latitude corresponds
    /// to the Y-coordinate remapped to a [0, 1] range from the north pole to the south pole.
    /// </summary>
    /// <param name="vertex">The 3D vertex on the unit sphere to convert to cylindrical UV coordinates.</param>
    /// <returns>A Vector2 representing the cylindrical UV coordinates, with longitude [0, 1] and latitude [0, 1].</returns>
    public static Vector2 CylindricalUV(this Vector3 vertex)
    {
        vertex.Normalize();

        // Longitude is the rotation about the north pole/y-axis (0-1)
        float longitude = Mathf.Atan2(vertex.z, vertex.x) * radToRots;
        // Remap -.5->+.5 to 0,1
        longitude += 0.5f;

        // Latitude is degrees from the north pole, towards the south pole (0-1)
        float latitude = vertex.y * 0.5f + 0.5f;

        return new Vector2(longitude, latitude);
    }

    /// <summary>
    /// Projects a 3D point onto a 2D plane defined by a normal and an origin point.
    /// The resulting 2D point is represented in the plane's local coordinate system.
    /// </summary>
    /// <param name="point">The 3D point to project onto the plane.</param>
    /// <param name="planeNormal">The normal vector defining the plane.</param>
    /// <param name="planeOrigin">A point on the plane used to define its position.</param>
    /// <returns>A 2D vector representing the projected point in the plane's local coordinate system.</returns>
    public static Vector2 ProjectPointOntoPlane(this Vector3 point, Vector3 planeNormal, Vector3 planeOrigin)
    {
        // Ensure planeNormal is not a zero vector
        if (planeNormal == Vector3.zero) throw new System.ArgumentException("ProjectPointOntoPlane failed: Plane normal cannot be zero.");

        Vector3 planeXAxis;
        Vector3 planeYAxis;

        if (planeNormal.CloseEqual(Vector3.up))
        {
            planeXAxis = Vector3.right;
            planeYAxis = Vector3.forward;  // Use a standard axis when normal is up
        }
        else
        {
            // Compute both axes with one cross product
            planeXAxis = Vector3.Cross(Vector3.up, planeNormal).normalized;
            planeYAxis = Vector3.Cross(planeNormal, planeXAxis);  // Only this second cross is needed
        }

        Vector3 offset = point - planeOrigin;
        // Project the point onto the plane
        Vector2 planePoint = new Vector2(
            Vector3.Dot(offset, planeXAxis),
            Vector3.Dot(offset, planeYAxis));

        return planePoint;
    }

    public class Vector3CloseComparer : IEqualityComparer<Vector3>
    {
        private readonly float tolerance;

        public Vector3CloseComparer(float tolerance = 0.0001f)
        {
            this.tolerance = tolerance;
        }

        public bool Equals(Vector3 v1, Vector3 v2)
        {
            return v1.CloseEqual(v2, tolerance);
        }

        /// <summary>
        /// Computes a hash code for a Vector3 based on a tolerance value for approximate equality.
        /// This method scales each component by the tolerance to create hash "buckets" for similar values.
        /// 
        /// Note:
        /// - Different Vector3 values may produce the same hash (collisions), which reduces efficiency.
        /// - Vector3 values considered equal by CloseEqual must produce the same hash for correctness.
        /// </summary>
        /// <param name="v">The Vector3 instance for which to compute the hash code.</param>
        /// <returns>An integer hash code that reflects the approximate equality of the Vector3.</returns>
        public int GetHashCode(Vector3 v)
        {
            // Scale components down by tolerance to bucket similar values
            int hashX = Mathf.RoundToInt(v.x / tolerance);
            int hashY = Mathf.RoundToInt(v.y / tolerance);
            int hashZ = Mathf.RoundToInt(v.z / tolerance);

            // Combine component hashes (e.g., using a tuple hash pattern)
            int hash = hashX;
            hash = (hash * 397) ^ hashY;
            hash = (hash * 397) ^ hashZ;
            return hash;
        }

    }

}

/// <summary>
/// Class used by the Polyhedron class to store the corners. Contains the actual vector3 position of the corner, and a reference to the polyhedron is is part of.
/// </summary>
[System.Serializable]
public class Corner 
{

    [HideInInspector]
    [System.NonSerialized]
    Polyhedron cornerOf;
    public void SetPoly(Polyhedron p)
    {
        cornerOf = p;
    }
    public Vector3 vertex;

    public Corner(Polyhedron cornerOf, Vector3 vertex)
    {
        this.cornerOf = cornerOf;
        this.vertex = vertex;
    }

    /// <summary>
    /// Potentially slow: looks though the polygon to find all the edges that touch this corner. Not cached.
    /// </summary>
    public List<Edge> touchingEdges
    {
        get 
        {
            List<Edge> edgesFound = new List<Edge>();
            foreach (Edge edgeToCheck in cornerOf.edges)
            {
                if (edgeToCheck.TouchesCorner(this))
                    edgesFound.Add(edgeToCheck);
            }
            Vector3 baseVector = edgesFound[0].VectorFrom(this);
            edgesFound.Sort((a, b) => 
                { 
                    if(Vector3.SignedAngle(baseVector, a.VectorFrom(this),vertex) >=
                        Vector3.SignedAngle(baseVector, b.VectorFrom(this), vertex))
                        return 1;
                    return -1;
                });
            return edgesFound;
        }
    }
    /// <summary>
    /// Potentially slow: looks though the polygon to find all the faces that touch this corner. Not cached.
    /// </summary>
    public List<Face> touchingFaces
    {
        get
        {
            List<Face> facesFound = new List<Face>();
            foreach (Face faceToCheck in cornerOf.faces)
            {
                if (faceToCheck.TouchesCorner(this))
                    facesFound.Add(faceToCheck);
            }
            /*facesFound.Sort((a, b) =>
            {
                if (Vector3.SignedAngle(facesFound[0].normal, a.normal, vertex) >
                    Vector3.SignedAngle(facesFound[0].normal, b.normal, vertex))
                    return 1;
                return -1;
            });*/
            return facesFound;
        }
    }
    public override string ToString()
    {
        return vertex.ToString();
    }
    /// <summary>
    /// Adds a new corner to the provided newCornersPloyhedron, for each corner in the provided cornersToCopy list.  This does NOT ensure that these corners are actually used by and lines or faces in the altered polyhedron.
    /// </summary>
    /// <param name="newCornersPloyhedron"></param>
    /// <param name="cornersToCopy"></param>
    /// <returns></returns>
    public static List<Corner> CopyCorners(Polyhedron newCornersPloyhedron, List<Corner> cornersToCopy)
    {
        List<Corner> newCornerList = new List<Corner>();
        foreach (Corner origCorner in cornersToCopy)
            newCornerList.Add(new Corner(newCornersPloyhedron, origCorner.vertex));
        return newCornerList;
    }

    /// <summary>
    /// Generates a list of Vector3's using the vectors referenced by each corner inthe provided corner list.
    /// </summary>
    /// <param name="cornerList">The list of corners from which to generate the vector3 list</param>
    /// <returns>the generated list of Vetcor3s</returns>
    public static List<Vector3> GetVerticies(List<Corner> cornerList)
    {
        List<Vector3> verts = new List<Vector3>();
        foreach (Corner c in cornerList)
            verts.Add(c.vertex);
        return verts;
    }


}

/// <summary>
/// a faces is made up of a list of corners
/// </summary>
[System.Serializable]
public class Face
{
    [HideInInspector]
    [System.NonSerialized]
    Polyhedron faceOf;
    public void SetPoly(Polyhedron p)
    {
        faceOf = p;
    }
    public List<Corner> corners;

    /// <summary>
    /// Creates a new face that references the provided polyhedron.  Does NOT add the face to the polyhedron.
    /// </summary>
    /// <param name="faceOf"></param>
    /// <param name="corners"></param>
    public Face(Polyhedron faceOf, List<Corner> corners)
    {
        this.faceOf = faceOf;
        this.corners = corners;
        if (corners.Count < 3)
            throw new GeometryException("Attempting to create a face with less that 3 corners: [" + string.Join(",", corners) + "]");
        for (int i = 0; i < 3; i++)
            if (corners[i] == null) throw new GeometryException("Attempting to create a face with null reference in corner: [" + i + "]");

        //confirm all corners lie on the same plane
        Vector3 planeNormal = normal;// Vector3.Cross(corners[1].vertex - corners[0].vertex, corners[2].vertex - corners[0].vertex).normalized;
        Vector3 planePoint = corners[0].vertex;
        for (int i = 3; i < corners.Count; i++) //triangles will skip this
        {
            //Vector3 testNormal = Vector3.Cross(corners[1].vertex - corners[0].vertex, corners[i].vertex - corners[0].vertex).normalized;
            Vector3 cornerToCheck = corners[i].vertex;
            float planeEq = //if point lies on the same plane- the result should be 0
                  (cornerToCheck.x - planePoint.x) * planeNormal.x
                + (cornerToCheck.y - planePoint.y) * planeNormal.y
                + (cornerToCheck.z - planePoint.z) * planeNormal.z;
            if (Mathf.Abs(planeEq) > .1f)
                Debug.Log("Attempting to add corner to a face, but corner does not lie on the same plane as the first three face corners. offset: " + planeEq);
            //throw new GeometryException("Attempting to add corner to a face, but corner does not lie on the same plane as the first three face corners. offset: " + planeEq);
            
        }
    }

    public List<Face> neighbors
    {
        get
        {
            List<Face> foundNeighbors = new List<Face>();
            foreach (Edge e in edges)
                foreach (Face f in e.touchingFaces)
                    if (f != this && !foundNeighbors.Contains(f))
                        foundNeighbors.Add(f);

            /*
            foreach (Face polyFace in faceOf.faces)
                if (polyFace != this)
                    foreach (Corner faceCorner in corners)
                        if (polyFace.TouchesCorner(faceCorner))
                        {
                            foundNeighbors.Add(polyFace);
                            break;// next face
                        }
            */
            return foundNeighbors;
        }
    }

    public Vector3 center
    {
        get
        {
            return Corner.GetVerticies(corners).AvgPos();
        }
    }

    public Vector3 normal
    {
        get { return Vector3.Cross(corners[1].vertex - corners[0].vertex, corners[2].vertex - corners[0].vertex).normalized; }
    }
        
    /// <summary>
    /// to do: cashe this- recompute on "corners" changed & in constructor/first pass
    /// </summary>
    public List<Edge> edges
    {
        get
        {
            List<Edge> edgesFound = new List<Edge>(2);

            foreach (Edge aPolyEdge in faceOf.edges)
            {
                int touchCount = 0;
                foreach (Corner aFaceCorner in corners)
                {
                    if (aPolyEdge.TouchesCorner(aFaceCorner))
                        touchCount++;
                }
                if (touchCount == 2)
                    edgesFound.Add(aPolyEdge);
            }

            return edgesFound;
        }
    }

    public bool TouchesEdge(Edge toCheck)
    {
        int sharedCornerCount = 0;
        foreach (Corner aFaceCorner in corners)
        {
            if (toCheck.TouchesCorner(aFaceCorner))
                sharedCornerCount++;
        }
        return sharedCornerCount == 2;
        
    }

    public bool TouchesCorner(Corner toCheck)
    {
        foreach (Corner aFaceCorner in corners)
            if (aFaceCorner == toCheck) return true;
        return false;        
    }


    public void ReOrderCornersClockWiseAroundCenterAndNormal()
    {
        //List<Corner> originalCorners = new List<Corner>(corners);// we will change the order of this member
        List<Vector3> faceVerts = Corner.GetVerticies(corners);
        Vector3 centerPt = faceVerts.AvgPos();
        Vector3 normal = centerPt.normalized;// this.normal;
        Vector3 startDir = corners[0].vertex - centerPt;
        corners.Sort((Corner a, Corner b) => 
            {
                Vector3 centerToA = a.vertex - centerPt;
                Vector3 centerToB = b.vertex - centerPt;
                if (Vector3.SignedAngle(startDir, centerToA,normal) < Vector3.SignedAngle(startDir, centerToB, normal))
                    return -1;
                return 1;
//                if (Vector3.Angle(startDir, centerToA) == Vector3.Angle(startDir, centerToB))
  //                  return 0;
            }
        );


    }

    public override string ToString()
    {
        string s=" Corners["+ corners.Count+"]:";
        foreach (Corner c in corners)
            s+= faceOf.corners.FindIndex((p) => { return c == p; }).ToString()+",";
        List<Edge> tempEdges= edges;
        s += "\n Edges["+ tempEdges .Count+ "]: ";
        foreach (Edge e in tempEdges)
            s += faceOf.edges.FindIndex((p) => { return e == p; }).ToString() + ",";
        return s;
    }
}
/// <summary>
/// An Edge specifies the border between two faces, and connects two corners.
/// </summary>
[System.Serializable]
public class Edge
{
    [HideInInspector]
    [System.NonSerialized]
    Polyhedron edgeOf;
    public void SetPoly(Polyhedron p)
    {
        edgeOf = p;
    }
    public Corner endpointA;
    public Corner endpointB;


    public Edge(Polyhedron edgeOf, Corner endpointA, Corner endpointB)
    {
        this.edgeOf = edgeOf;
        this.endpointA = endpointA;
        this.endpointB = endpointB;
    }

    public float length
    {
        get { return (endpointB.vertex - endpointA.vertex).magnitude; }
    }
    // returns false is the edge no longer exists
    public bool ReplaceCorner(Corner cornerToReplace, Corner replacement)
    {
        if (endpointA == cornerToReplace)
            endpointA = replacement;
        else if (endpointB == cornerToReplace)
            endpointB = replacement;
        else
            throw new GeometryException("Trying to replace a corner on an Edge, but the edge does not contain the corner specified.");
        return (endpointA != endpointB);
    }


    public bool TouchesCorner(Corner toCheck)
    {
        return toCheck == endpointA || toCheck == endpointB;
    }

    public List<Face> touchingFaces
    {
        get {
            List<Face> touchingFaces = new List<Face>();
            foreach (Face existingFace in edgeOf.faces) // find both faces that touch the edge being truncated
            {
                if (existingFace.TouchesEdge(this))
                    touchingFaces.Add(existingFace);
            }
            if (touchingFaces.Count > 2)
                throw new GeometryException("Finding faces touching edge- but more that 2 faces found.");
            if (touchingFaces.Count < 2)
                Debug.Log("Warning: Finding faces touching edge, but it appears to be touching less than 2 faces.  2d-objects can ignore this warning.");
            return touchingFaces;
        }
    }

    public Vector3 VectorFrom(Corner fromCorner)
    {
        if (fromCorner == endpointA)
        {
            return (endpointB.vertex - endpointA.vertex);
        }
        if (fromCorner == endpointB)
        {
            return (endpointA.vertex - endpointB.vertex);
        }
        throw new GeometryException("Trying to find vector from corner, on an edge; but corner passed to function is NOT on the edge.");
    }

    public Vector3 FractionFrom(float fraction, Corner fromCorner)
    {
        if (fromCorner == endpointA)
        {
            return (endpointB.vertex - endpointA.vertex) * fraction + endpointA.vertex;
        }
        if (fromCorner == endpointB)
        {
            return (endpointA.vertex - endpointB.vertex) * fraction + endpointB.vertex;
        }
        throw new GeometryException("Trying to find fraction from corner, on an edge; but corner passed to function is NOT on the edge.");
    }


    public Vector3 midPoint
    {
        get {
            return (endpointA.vertex + endpointB.vertex) * 0.5f;
        }
    }

    // Override Equals to define equality based on endpoints
    public override bool Equals(object obj)
    {
        if (obj is Edge otherEdge)
        {
            return (endpointA == otherEdge.endpointA && endpointB == otherEdge.endpointB) ||
                   (endpointA == otherEdge.endpointB && endpointB == otherEdge.endpointA);
        }
        return false;
    }

    // Override GetHashCode to ensure equal objects have the same hash code
    public override int GetHashCode()
    {
        int hash1 = endpointA.GetHashCode() ^ endpointB.GetHashCode();
        int hash2 = endpointB.GetHashCode() ^ endpointA.GetHashCode();
        return hash1 ^ hash2;
    }

    // Override == and != operators for usability
    public static bool operator ==(Edge e1, Edge e2)
    {
        if (ReferenceEquals(e1, e2)) return true;
        if (e1 is null || e2 is null) return false;
        return e1.Equals(e2);
    }

    public static bool operator !=(Edge e1, Edge e2)
    {
        return !(e1 == e2);
    }
    public override string ToString()
    {
        
        int indexA = edgeOf.corners.FindIndex((c) => { return c == endpointA; });
        int indexB = edgeOf.corners.FindIndex((c) => { return c == endpointB; });
        return "Edge- endpointA["+indexA+"]: " + endpointA + ", endpointB[" + indexB + "]: " + endpointB;
    }


}


[System.Serializable]
public class Polyhedron : ISerializationCallbackReceiver
{
    
    public List<Face> faces;
    public List<Edge> edges;
    public List<Corner> corners; // only corners store the Vector3 positions.  Corners are shared between faces and between edges.
    
    public Polyhedron(Polyhedron toCopy)
    {
        this.corners = new List<Corner>(toCopy.corners);
        this.edges = new List<Edge>(toCopy.edges);
        this.faces = new List<Face>(toCopy.faces);


        foreach (Corner c in corners)
            c.SetPoly(this);
        foreach (Edge c in edges)
            c.SetPoly(this);
        foreach (Face c in faces)
            c.SetPoly(this);
    }

    public Polyhedron(List<Vector3> cornerPos, List<List<int>> facesByCornerIndex)
    {
        corners = new List<Corner>();
        foreach (Vector3 pos in cornerPos)
            corners.Add(new Corner(this, pos));
        faces = new List<Face>();
        edges = new List<Edge>();
        foreach (List<int> faceCornerByIndexList in facesByCornerIndex)
        {
            //generate corner list from index lest
            List<Corner> faceCorners = new List<Corner>();
            foreach (int index in faceCornerByIndexList)
                faceCorners.Add(corners[index]);
            //create face
            faces.Add(new Face(this, faceCorners));

        }// end loop all faces
        RecomputeEdges();
    }

    public Polyhedron(Mesh fromMesh)
    {
        corners = new List<Corner>(); //use  new public Corner(Polyhedron cornerOf, Vector3 vertex) for each
        faces = new List<Face>();  //use new Face(Polyhedron faceOf, List<Corner> corners)
        edges = new List<Edge>(); //use new Edge(Polyhedron edgeOf, Corner endpointA, Corner endpointB)

        //loop through mesh triangles and vertex array to populate the coner face and edge lists- start with corners
        Dictionary<Vector3, Corner> cornerByPosition = new Dictionary<Vector3, Corner>();
        int[] vertIndexToCornerIndex = new int[fromMesh.vertexCount]; // provide mapping from vertex index into a corner index (in the corners list- only valid until corners are removed from this list, later on)
        for (int vertIndex = 0; vertIndex < fromMesh.vertexCount; vertIndex++)
        {
            Vector3 vertPos = fromMesh.vertices[vertIndex];
            if (!cornerByPosition.ContainsKey(vertPos))
            {
                Corner newCorner = new Corner(this, fromMesh.vertices[vertIndex]);
                corners.Add(newCorner);
                cornerByPosition.Add(vertPos, newCorner);
                vertIndexToCornerIndex[vertIndex] = corners.Count-1;
            }// skip if vert is already a corner
            else
                vertIndexToCornerIndex[vertIndex] = corners.FindIndex((c)=> { return c == cornerByPosition[vertPos]; });
        }

        //we will start with one face per triangle, then combine flush faces that are on the same plane- eliminating the edges between them.
        // Helper function to find or create an edge between two corners
        Edge GetOrCreateEdge(Corner endpointA, Corner endpointB)
        {
            // Check if an edge already exists between these corners
            foreach (var edge in edges)
            {
                if ((edge.endpointA == endpointA && edge.endpointB == endpointB) ||
                    (edge.endpointA == endpointB && edge.endpointB == endpointA))
                {
                    return edge;
                }
            }

            // If not, create a new edge and add it to the edges list
            Debug.Log("new edge created ["+edges.Count+"] A:"+endpointA+"  B:"+ endpointB);
            var newEdge = new Edge(this, endpointA, endpointB);
            edges.Add(newEdge);
            return newEdge;
        }

        // Step 2: Loop through each triangle in the mesh to create faces
        for (int triIndex = 0; triIndex < fromMesh.triangles.Length; triIndex += 3)
        {
            // Get the three vertex indices of the triangle
            int vertAIndex = fromMesh.triangles[triIndex];
            int vertBIndex = fromMesh.triangles[triIndex + 1];
            int vertCIndex = fromMesh.triangles[triIndex + 2];

            // Retrieve or create the corners for this triangle
            Corner cornerA = corners[vertIndexToCornerIndex[vertAIndex]];
            Corner cornerB = corners[vertIndexToCornerIndex[vertBIndex]];
            Corner cornerC = corners[vertIndexToCornerIndex[vertCIndex]];

            // Create the face with its corners and add it to the face list
            List<Corner> faceCorners = new List<Corner> { cornerA, cornerB, cornerC };
            Face face = new Face(this, faceCorners);
            
            faces.Add(face);
        }
        RecomputeEdges();
        bool FacesShareAnEdge(List<Edge> faceAedges, List<Edge>faceBedges, out Edge sharedEdgeA, out Edge sharedEdgeB)
        {
          
            foreach (Edge edgeA in faceAedges)
            {
                foreach (Edge edgeB in faceBedges)
                {
                    if (edgeA == edgeB)
                    {
                        sharedEdgeA = edgeA;
                        sharedEdgeB = edgeB;
                     //   Debug.Log("..checking forsharedEdges found: " + sharedEdgeA);
                        return true;
                    }
                }
            }
            sharedEdgeA = null;
            sharedEdgeB = null;
           // Debug.Log("..checking forsharedEdges found: none");
            return false;
        }
        bool FacesArePlanar(Face faceA, Face faceB)
        {
          //  Debug.Log("..checking normals: " + faceA.normal + " , " + faceB.normal +" :  "+ Vector3.Dot(faceA.normal, faceB.normal));
            return Mathf.Abs(Vector3.Dot(faceA.normal,faceB.normal)) > .99f;
        }
        //merge faces
        bool noMergesDone = false;
        while (!noMergesDone)
        {
            noMergesDone = true;
            for (int faceIndex = 0; faceIndex < faces.Count; faceIndex++)
            {
                Face faceA = faces[faceIndex];
                List<Edge> faceAEdges = faceA.edges; // this is the faceobject we will keep,the other will be discarded
                for (int otherFaceIndex = faceIndex + 1; otherFaceIndex < faces.Count; otherFaceIndex++)
                {
                    Face faceB = faces[otherFaceIndex];
                    List<Edge> faceBEdges = faceB.edges;
                    // Debug.Log("checking faces for merge: " + faceIndex + " , " + otherFaceIndex);
                    if (FacesArePlanar(faceA, faceB) && FacesShareAnEdge(faceA.edges, faceB.edges, out Edge ignore, out Edge ignoreB))
                    {
                        //Debug.Log("Merging faces: " + faceIndex + " , " + otherFaceIndex + "\n" + faceA + "\n" + faceB);
                        
                        noMergesDone = false;
                        //combine faces
                       
                        
                        //remove all shared edges from BOTH faces
                        while (FacesShareAnEdge(faceAEdges, faceBEdges, out Edge sharedEdgeA, out Edge sharedEdgeB))
                        {
                            //Debug.Log("  shared edge["+ faceIndex + "]: " + sharedEdgeA + ", and edge[" + faceIndex + "]: " + sharedEdgeB);
                            faceAEdges.Remove(sharedEdgeA);//remove edge from face
                            faceBEdges.Remove(sharedEdgeB);//remove edge from face
                                                           // Debug.Log("Removing shared edge: " + sharedEdge + " edges.count: " + edges.Count);

                           // edges.Remove(sharedEdgeA);//remove edge from polyhedron
                           // if (!ReferenceEquals(sharedEdgeA, sharedEdgeB))
                             //   edges.Remove(sharedEdgeB);//remove edge from polyhedron
                          //  Debug.Log("Removed shared edge: " + sharedEdgeA + " edges.count: " + edges.Count);

                        }
                        //add remaining edges in B to A
                        foreach (Edge e in faceBEdges)
                        {
                            if (!faceAEdges.Contains(e))
                                faceAEdges.Add(e);
                        }

                        //recompute faceA corners from all edges it now has
                        faceA.corners = new List<Corner>();
                        foreach (Edge e in faceAEdges)
                        {
                            if (!faceA.corners.Contains(e.endpointA)) faceA.corners.Add(e.endpointA);
                            if (!faceA.corners.Contains(e.endpointB)) faceA.corners.Add(e.endpointB);
                        }
                        //Debug.Log("  face ["+faceIndex+"]  now has " + faceA.corners.Count + " corners and " + faceAEdges.Count + " edges");

                        faces.RemoveAt(otherFaceIndex);// remove face from polyhedron
                       // Debug.Log("Removed face at index: " + otherFaceIndex + " faces.count: " + faces.Count);
                        otherFaceIndex--;// undo ++ in for loop

                    }//end faces are planar
                    else
                    {
                       // Debug.Log("faces not flush and parallel");
                    }
                }//loop otherfaces
                faceA.ReOrderCornersClockWiseAroundCenterAndNormal();
            }//loop all faces
            //Debug.Log("MergePassDone. any merges:"+ !noMergesDone);
        }// all merges done
        //recompute poly corner listto onlyinclude those in remainig faces
        
        List<Corner> newCorners = new List<Corner>();
        //List<Edge> newEdges = new List<Edge>();

       // Debug.Log("Rebuilding from faces: " + faces.Count);
        int faceCount = 0;
        foreach (Face face in faces)
        {
            //Debug.Log(" +face ["+faceCount+"] corners: " + (face.corners.Count) );
            faceCount++;
            foreach (Corner corner in face.corners)
            {
                if(!newCorners.Contains(corner))
                    newCorners.Add(corner);
            }
           // Debug.Log(" =total(no dup) corners,edges: " + newCorners.Count + "," + newEdges.Count);
        }
        corners = newCorners;
        RecomputeEdges();
        //edges = newEdges;
    }// end contructor from mesh
    //empty
    public Polyhedron()
    {
        corners = new List<Corner>();
        faces = new List<Face>();
        edges = new List<Edge>();
    }

    public class CancelBoolRef
    {
        public bool doCancel = false;
    }
    class YieldTimer
    {
        private readonly int timeSlice; // Time slice in milliseconds
        private readonly System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
        private readonly CancelBoolRef cancelRef; // Optional cancellation reference

        public YieldTimer(int timeSlice = 10, CancelBoolRef cancelRef = null)
        {
            this.timeSlice = timeSlice;
            this.cancelRef = cancelRef;
            timer.Start();
        }
        public YieldTimer(CancelBoolRef cancelRef)
        {
            this.timeSlice = 10;
            this.cancelRef = cancelRef;
            timer.Start();
        }
        public async UniTask YieldOnTimeSlice()
        {
            if (cancelRef != null && cancelRef.doCancel)
                throw new System.OperationCanceledException();
            if (timer.ElapsedMilliseconds > timeSlice)
            {
                await UniTask.Yield(PlayerLoopTiming.Update);
                timer.Restart(); // Restart timer after yielding
            }
        }
    }

    /// <summary>
    /// After creating the polygon by specifying corners and faces, the function may be called to recompute all edges based upon this information.
    /// Ensures only one edge exists between each pair of touching faces.
    /// </summary>
    async UniTask RecomputeEdgesAsync(CancelBoolRef cancelRef)
    {
        YieldTimer yieldTimer = new YieldTimer(cancelRef);
        edges.Clear();
        foreach (Face eachFace in faces)
        {

            Vector3 faceNormal = eachFace.normal;
            Vector3 faceCenter = eachFace.center;
            Vector3 vectorToFirstCorner = eachFace.corners[0].vertex - faceCenter;
            eachFace.corners.Sort((a, b) =>
            {
                if (Vector3.SignedAngle(vectorToFirstCorner, a.vertex - faceCenter, faceNormal)
                    >= Vector3.SignedAngle(vectorToFirstCorner, b.vertex - faceCenter, faceNormal))
                    return 1;
                return -1;
            });
            await yieldTimer.YieldOnTimeSlice();
            /*
            //debug stuff
            Debug.Log("JustSorted "+eachFace.corners.Count+"face Corners  face center "+ faceCenter);
            foreach (Corner c in eachFace.corners)
            {
                float angle = Vector3.SignedAngle(vectorToFirstCorner, c.vertex - faceCenter, faceNormal);
                Debug.Log("    Sorted corner on face at pos from center "+ (c.vertex - faceCenter) + ", angle from first corner: " + angle);
            }

            // end debug stuff
            */
            for (int i = 1; i <= eachFace.corners.Count; i++)
            {
                Corner endpointA = eachFace.corners[i - 1];//start point:  i-1 = zero on first pass
                Corner endpointB;
                if (i == eachFace.corners.Count)//be circular
                    endpointB = eachFace.corners[0];//endpoint: will be zero on last pass
                else
                    endpointB = eachFace.corners[i];//  will be 1 through eachFace.corners.Count-1 (on every pass but last)


                bool doAdd = true;
                foreach (Edge testEdge in edges)
                    if (testEdge.TouchesCorner(endpointA) && testEdge.TouchesCorner(endpointB))// it touches BOTH corners- same edge
                    {
                        doAdd = false;
                        break;// no need to keep looping existing edges
                    }
                if (doAdd)
                    edges.Add(new Edge(this, endpointA, endpointB));
                await yieldTimer.YieldOnTimeSlice();
            }
        }

    }
    void RecomputeEdges()
    {
        edges.Clear();
        foreach (Face eachFace in faces)
        {
            
            Vector3 faceNormal = eachFace.normal;
            Vector3 faceCenter = eachFace.center;
            Vector3 vectorToFirstCorner = eachFace.corners[0].vertex - faceCenter;
            eachFace.corners.Sort((a, b) => 
                {
                    if (Vector3.SignedAngle(vectorToFirstCorner, a.vertex- faceCenter, faceNormal)
                        >= Vector3.SignedAngle(vectorToFirstCorner, b.vertex- faceCenter, faceNormal))
                        return 1;
                    return -1; 
                });
            /*
            //debug stuff
            Debug.Log("JustSorted "+eachFace.corners.Count+"face Corners  face center "+ faceCenter);
            foreach (Corner c in eachFace.corners)
            {
                float angle = Vector3.SignedAngle(vectorToFirstCorner, c.vertex - faceCenter, faceNormal);
                Debug.Log("    Sorted corner on face at pos from center "+ (c.vertex - faceCenter) + ", angle from first corner: " + angle);
            }

            // end debug stuff
            */
            for (int i = 1; i <= eachFace.corners.Count; i++)
            {
                Corner endpointA = eachFace.corners[i - 1];//start point:  i-1 = zero on first pass
                Corner endpointB;
                if (i == eachFace.corners.Count)//be circular
                    endpointB = eachFace.corners[0];//endpoint: will be zero on last pass
                else
                    endpointB = eachFace.corners[i];//  will be 1 through eachFace.corners.Count-1 (on every pass but last)


                bool doAdd = true;
                foreach (Edge testEdge in edges)
                    if (testEdge.TouchesCorner(endpointA) && testEdge.TouchesCorner(endpointB))// it touches BOTH corners- same edge
                    {
                        doAdd = false;
                        break;// no need to keep looping existing edges
                    }
                if (doAdd)
                    edges.Add(new Edge(this, endpointA, endpointB));
            }
        }

    }

    Corner CreateOrFindCornerByVertex(Vector3 vertex)
    {
        foreach (Corner c in corners)
            if ((vertex-c.vertex).sqrMagnitude<.01f)
                return c;
        return new Corner(this, vertex);
    }


    Dictionary<Vector3, Corner> cornerByVectorCache = new Dictionary<Vector3, Corner>(new Vector3Extensions.Vector3CloseComparer());
    /// <summary>
    /// Initializes the cache for mapping vertices to corners. 
    /// Ensures no duplicate vertices exist in the polyhedron.
    /// </summary>
    /// <exception cref="GeometryException">
    /// Thrown when duplicate vertices are detected in the polyhedron.
    /// </exception>
    void InitCornerByVectorCache()
    {
        cornerByVectorCache.Clear();
        foreach (Corner c in corners)
        {
            if (cornerByVectorCache.TryGetValue(c.vertex, out Corner corner))
            {
                throw new GeometryException("Unexpected duplicate Vector value found in polyhedron:" + c.vertex);
            }
            cornerByVectorCache.Add(c.vertex, c);
        }
    }
    /// <summary>
    /// Searches all corners in the polyhedron for one that specifies the provided vertex.
    /// Utilizes a cached dictionary for efficient lookup; falls back to a linear search if necessary.  Success will add the found corner, and it's vertex
    /// to the cache  
    /// </summary>
    /// <param name="vertex">The vertex to locate.</param>
    /// <returns>
    /// The matching Corner object, or null if no match is found.
    /// </returns>
    Corner FindCornerByVertex(Vector3 vertex)
    {
        if (cornerByVectorCache.TryGetValue(vertex, out Corner corner))
            return corner;
        foreach (Corner c in corners)
        {
            if (vertex.CloseEqual(c.vertex))
            {
                cornerByVectorCache.Add(c.vertex, c);
                return c;
            }
        }
        return null;
    }


    float AverageEdgeLength(Corner c)
    {
        float totalLength = 0f;
        int edgeCount = c.touchingEdges.Count;
        for (int i = 0; i < edgeCount; i++)
        {
            totalLength += c.touchingEdges[i].length;
        }
        return edgeCount > 0 ? totalLength / edgeCount : 0f;
    }

    float UniformEdgeLengthTruncFraction(Corner c)//, Edge edgeToTruncate)
    {

        Vector3 firstEdge = c.touchingEdges[0].VectorFrom(c);
        Vector3 secondEdge = c.touchingEdges[1].VectorFrom(c);
        float edgeLen = firstEdge.magnitude;
        float angleBetween = Mathf.Deg2Rad * Vector3.Angle(firstEdge, secondEdge);


        float halfOppLen = Mathf.Sin(.5f * angleBetween);
        return halfOppLen / (halfOppLen + 1);

    }

    float UniformFaceCenterToConerLengthTruncFraction(Corner c)
    {
        Vector3 firstEdge = c.touchingEdges[0].VectorFrom(c);
        Vector3 secondEdge = c.touchingEdges[1].VectorFrom(c);
        float edgeLen = firstEdge.magnitude;
        float angleBetween = Mathf.Deg2Rad * Vector3.Angle(firstEdge, secondEdge);
        Vector3 oldFaceCenter = c.touchingEdges[0].touchingFaces[0].center;

        Vector3 oldCornerToFaceCenter = c.vertex - oldFaceCenter;

        float M= oldCornerToFaceCenter.magnitude;
        float m = angleBetween * .5f;

        float C = Mathf.Sin(m) * M;

        float angleD = Vector3.Angle(-c.vertex, firstEdge) * Mathf.Deg2Rad;
        float d = angleD;
        float sd = Mathf.Sin(d);
        float K = edgeLen * 0.5f;

        float secd = 1f / Mathf.Cos(d);
        float B = -(secd * secd) * (Mathf.Sqrt((C * C + K * K) * sd * sd - C * C) - K);

        return B / edgeLen;

    }

    float UniformFaceCenterToEdgeMidpoint(Corner c)
    {
        Edge firstEdgeRef = c.touchingEdges[0];
        Vector3 firstEdge = firstEdgeRef.VectorFrom(c);
        Vector3 secondEdge = c.touchingEdges[1].VectorFrom(c);

        Vector3 oldFaceCenter = firstEdgeRef.touchingFaces[0].center;
        Vector3 oldFaceCenterToEdgeMid = firstEdgeRef.midPoint - oldFaceCenter;
        Vector3 oldFaceCenterToCorner = c.vertex-oldFaceCenter;

        float angleC = Vector3.Angle(-c.vertex, oldFaceCenterToCorner) * Mathf.Deg2Rad;
        float angleR = Vector3.Angle(firstEdge, secondEdge) * Mathf.Deg2Rad * 0.5f;
        //B = C/(sin(c)*cos(r))
        float B = oldFaceCenterToEdgeMid.magnitude / (firstEdge.magnitude*Mathf.Sin(angleC) * Mathf.Cos(angleR));
        return B;

    }


    //cuts the polyhedron such that each corner becomes a face (removing the old corner, and replacing it with new corners)
    //the number of corners this new face will have is equal to the number of edges touching the original corner.
    //existing faces that touched the original corner will end up with TWO new corners, and a new edge, and will NOT have the original corner anymore.
    public async UniTask<Polyhedron> TruncAsync(float depthAsFractionOfEdgeLength = 0.5f,CancelBoolRef cancelRef=null)
    {
        YieldTimer timer = new YieldTimer(cancelRef);
        Polyhedron truncPoly = new Polyhedron();
        List<Corner> newCorners = truncPoly.corners;//  this list- after first loop- we be the same size as- and in the same order as ``this.faces``
        List<Face> newFaces = truncPoly.faces;

        // The index of these corners lists will match those of faces in original poly
        // we will populate each of this as we loop through old corner's edges
        // after the loop we will use them to create the new versions of the old faces.
        List<List<Corner>> newReplacedFaceCorners = new List<List<Corner>>();
        foreach (Face oldFace in this.faces)
        {
            newReplacedFaceCorners.Add(new List<Corner>());
        }

        foreach (Corner oldCorner in this.corners)
        {
            depthAsFractionOfEdgeLength = UniformEdgeLengthTruncFraction(oldCorner);

            //corner being chopped of- this will be the corners of the face that remains
            List<Corner> newFaceCorners = new List<Corner>();

            //    Debug.Log("looping old corner- touching "+ oldCorner.touchingEdges.Count + " edges");
            // all edges that touch this corner
            foreach (Edge oldEdgeTouchingOldCorner in oldCorner.touchingEdges)
            {
                //  Debug.Log("DepthFraction: " + depthAsFractionOfEdgeLength);
                Vector3 newCornerPos = oldEdgeTouchingOldCorner.FractionFrom(depthAsFractionOfEdgeLength, oldCorner);
                Corner newCorner = truncPoly.FindCornerByVertex(newCornerPos);
                if (newCorner == null)
                {
                    newCorner = new Corner(truncPoly, newCornerPos);
                    truncPoly.corners.Add(newCorner);
                }
                newFaceCorners.Add(newCorner);

                // old faces that touch this edge need to use the new corner instead of the old one now
                // so we add it to the list at newReplacedFaceCorners[oldFaceIndex]
                foreach (Face oldFaceTouchingCorner in oldEdgeTouchingOldCorner.touchingFaces)
                {
                    int oldFaceIndex = this.faces.IndexOf(oldFaceTouchingCorner);
                    newReplacedFaceCorners[oldFaceIndex].Add(newCorner);
                    //      Debug.Log("added corner to list at old face index: " + oldFaceIndex + "  List now contains " +
                    //          newReplacedFaceCorners[oldFaceIndex].Count + " corners");
                }
            }
            // Debug.Log("adding NEW face with " + newFaceCorners.Count+ " corners");
            truncPoly.faces.Add(new Face(truncPoly, newFaceCorners)); // new face to replace oldCorner
            await timer.YieldOnTimeSlice();
        }
        //add in original faces- "now with new corners!"
        foreach (List<Corner> newReplacementFaceCornerList in newReplacedFaceCorners)
        {
            truncPoly.faces.Add(new Face(truncPoly, newReplacementFaceCornerList));
            await timer.YieldOnTimeSlice();
        }

        foreach (Face f in newFaces)
        {
            f.ReOrderCornersClockWiseAroundCenterAndNormal();
            await timer.YieldOnTimeSlice();
        }

        await truncPoly.RecomputeEdgesAsync(cancelRef);
        // truncPoly.CheckIntegrity();
        return truncPoly;
    }
    public Polyhedron Trunc(float depthAsFractionOfEdgeLength = 0.5f)
    {
        Polyhedron truncPoly = new Polyhedron();
        List<Corner> newCorners = truncPoly.corners;//  this list- after first loop- we be the same size as- and in the same order as ``this.faces``
        List<Face> newFaces = truncPoly.faces;

        // The index of these corners lists will match those of faces in original poly
        // we will populate each of this as we loop through old corner's edges
        // after the loop we will use them to create the new versions of the old faces.
        List<List<Corner>> newReplacedFaceCorners = new List<List<Corner>>();
        foreach (Face oldFace in this.faces)
        {
            newReplacedFaceCorners.Add(new List<Corner>());
        }

        foreach (Corner oldCorner in this.corners)
        {
            depthAsFractionOfEdgeLength = UniformEdgeLengthTruncFraction(oldCorner);

            //corner being chopped of- this will be the corners of the face that remains
            List<Corner> newFaceCorners = new List<Corner>();

        //    Debug.Log("looping old corner- touching "+ oldCorner.touchingEdges.Count + " edges");
            // all edges that touch this corner
            foreach (Edge oldEdgeTouchingOldCorner in oldCorner.touchingEdges)
            {
              //  Debug.Log("DepthFraction: " + depthAsFractionOfEdgeLength);
                Vector3 newCornerPos = oldEdgeTouchingOldCorner.FractionFrom(depthAsFractionOfEdgeLength, oldCorner);
                Corner newCorner = truncPoly.FindCornerByVertex(newCornerPos);
                if (newCorner == null)
                {
                    newCorner = new Corner(truncPoly, newCornerPos);
                    truncPoly.corners.Add(newCorner);
                }
                newFaceCorners.Add(newCorner);

                // old faces that touch this edge need to use the new corner instead of the old one now
                // so we add it to the list at newReplacedFaceCorners[oldFaceIndex]
                foreach (Face oldFaceTouchingCorner in oldEdgeTouchingOldCorner.touchingFaces)
                {
                    int oldFaceIndex = this.faces.IndexOf(oldFaceTouchingCorner);
                    newReplacedFaceCorners[oldFaceIndex].Add(newCorner);
              //      Debug.Log("added corner to list at old face index: " + oldFaceIndex + "  List now contains " +
              //          newReplacedFaceCorners[oldFaceIndex].Count + " corners");
                }
            }
           // Debug.Log("adding NEW face with " + newFaceCorners.Count+ " corners");
            truncPoly.faces.Add(new Face(truncPoly, newFaceCorners)); // new face to replace oldCorner
        }
        //add in original faces- "now with new corners!"
        foreach (List<Corner> newReplacementFaceCornerList in newReplacedFaceCorners)
        {
            truncPoly.faces.Add(new Face(truncPoly, newReplacementFaceCornerList));
        }

        foreach (Face f in newFaces)
            f.ReOrderCornersClockWiseAroundCenterAndNormal();

        truncPoly.RecomputeEdges();
       // truncPoly.CheckIntegrity();
        return truncPoly;
    }

    public async UniTask<Polyhedron> TesselateFacesRadialAsync(CancelBoolRef cancelRef)
    {
        YieldTimer yieldTimer = new YieldTimer(cancelRef);
        Polyhedron tessPoly = new Polyhedron(this);
        List<Face> replacementFaces = new List<Face>();
        foreach (Face f in tessPoly.faces)
        {
            Corner newCenterCorner = new Corner(tessPoly, f.center);
            tessPoly.corners.Add(newCenterCorner);
            foreach (Edge e in f.edges) // each edge will become a new (triangular) face that touches centerpt
            {
                Face newFace = new Face(tessPoly, new List<Corner>() { e.endpointA, e.endpointB, newCenterCorner });
                newFace.ReOrderCornersClockWiseAroundCenterAndNormal();
                replacementFaces.Add(newFace);
                await yieldTimer.YieldOnTimeSlice();
            }
        }
        tessPoly.faces = replacementFaces;
        await tessPoly.RecomputeEdgesAsync(cancelRef);
        return tessPoly;
    }
    public Polyhedron TesselateFacesRadial()
    {
        Polyhedron tessPoly = new Polyhedron(this);
        List<Face> replacementFaces = new List<Face>();
        foreach (Face f in tessPoly.faces)
        {
            Corner newCenterCorner = new Corner(tessPoly,f.center);
            tessPoly.corners.Add(newCenterCorner);
            foreach (Edge e in f.edges) // each edge will become a new (triangular) face that touches centerpt
            {
                Face newFace = new Face(tessPoly, new List<Corner>() { e.endpointA, e.endpointB, newCenterCorner });
                newFace.ReOrderCornersClockWiseAroundCenterAndNormal();
                replacementFaces.Add(newFace);
            }
        }
        tessPoly.faces = replacementFaces;
        tessPoly.RecomputeEdges();
        return tessPoly;
    }

    public Polyhedron TesselateTriangleByEdgeMiddles()
    {
        Polyhedron tessPoly = new Polyhedron();
        InitCornerByVectorCache();

        foreach (Face f in faces)
        {
            if (f.corners.Count > 3)
                throw new GeometryException("The TesselateTriangleByEdgeMiddles, requires that all faces are triangles, but a face with " + f.corners.Count + " corners was found.");

            // corner ABC
            // mids   AB, BC, CA
            // new tris
            //  A,AB,AC
            //  B,BC,AB
            //  C,CA,BC
            //  AB,BC,CA
            Vector3 A = f.corners[0].vertex;
            Vector3 B = f.corners[1].vertex;
            Vector3 C = f.corners[2].vertex;

            Vector3 AB = (A + B) * 0.5f;
            Vector3 BC = (B + C) * 0.5f;
            Vector3 CA = (C + A) * 0.5f;

            //popout midpoints (no, we have a spherize function for that)
            /*float mag = A.magnitude;
            AB = AB.normalized * mag;
            BC = BC.normalized * mag;
            CA = CA.normalized * mag;
            */

            //get/create corners
            Corner cornerA = tessPoly.FindCornerByVertex(A);
            if (cornerA == null)
            {
                cornerA = new Corner(tessPoly, A);
                tessPoly.corners.Add(cornerA);
            }
            Corner cornerB = tessPoly.FindCornerByVertex(B);
            if (cornerB == null){
                cornerB = new Corner(tessPoly, B);
                tessPoly.corners.Add(cornerB);
            }
            Corner cornerC = tessPoly.FindCornerByVertex(C);
            if (cornerC == null)
            {
                cornerC = new Corner(tessPoly, C);
                tessPoly.corners.Add(cornerC);
            }

            Corner cornerAB = tessPoly.FindCornerByVertex(AB);
            if (cornerAB == null){
                cornerAB = new Corner(tessPoly, AB);
                tessPoly.corners.Add(cornerAB);
            }

            Corner cornerBC = tessPoly.FindCornerByVertex(BC);
            if (cornerBC == null)
            {
                cornerBC = new Corner(tessPoly, BC);
                tessPoly.corners.Add(cornerBC);
            }

            Corner cornerCA = tessPoly.FindCornerByVertex(CA);
            if (cornerCA == null)
            {
                cornerCA = new Corner(tessPoly, CA);
                tessPoly.corners.Add(cornerCA);
            }

            Face[] newFaces = new Face[4];
            newFaces[0] = new Face(tessPoly, new List<Corner>()
                {
                    cornerA,
                    cornerAB,
                    cornerCA
                });
            newFaces[1] = new Face(tessPoly, new List<Corner>()
                {
                    cornerB,
                    cornerBC,
                    cornerAB
                });
            newFaces[2] = new Face(tessPoly, new List<Corner>()
                {
                    cornerC,
                    cornerCA,
                    cornerBC
                });
            newFaces[3] = new Face(tessPoly, new List<Corner>()
                {
                    cornerAB,
                    cornerBC,
                    cornerCA
                });
            
            tessPoly.faces.AddRange(newFaces);

        }// end loop all faces

        tessPoly.RecomputeEdges();
        return tessPoly;
    }

    
    public async UniTask<Polyhedron> TesselateTriangleByEdgeMiddlesAsync(CancelBoolRef cancelRef)
    {
        Polyhedron tessPoly = new Polyhedron();
        InitCornerByVectorCache();

        YieldTimer yieldTimer = new YieldTimer(cancelRef);
        foreach (Face f in faces)
        {
            if (cancelRef.doCancel)
                throw new System.Exception("Process canceled");

            if (f.corners.Count > 3)
                throw new GeometryException("The TesselateTriangleByEdgeMiddlesAsync requires that all faces are triangles, but a face with " + f.corners.Count + " corners was found.");

            Vector3 A = f.corners[0].vertex;
            Vector3 B = f.corners[1].vertex;
            Vector3 C = f.corners[2].vertex;

            Vector3 AB = (A + B) * 0.5f;
            Vector3 BC = (B + C) * 0.5f;
            Vector3 CA = (C + A) * 0.5f;

            Corner cornerA = await FindOrCreateCornerAsync(tessPoly, A);
            Corner cornerB = await FindOrCreateCornerAsync(tessPoly, B);
            Corner cornerC = await FindOrCreateCornerAsync(tessPoly, C);

            Corner cornerAB = await FindOrCreateCornerAsync(tessPoly, AB);
            Corner cornerBC = await FindOrCreateCornerAsync(tessPoly, BC);
            Corner cornerCA = await FindOrCreateCornerAsync(tessPoly, CA);

            Face[] newFaces = new Face[4];
            newFaces[0] = new Face(tessPoly, new List<Corner>()
        {
            cornerA,
            cornerAB,
            cornerCA
        });
            newFaces[1] = new Face(tessPoly, new List<Corner>()
        {
            cornerB,
            cornerBC,
            cornerAB
        });
            newFaces[2] = new Face(tessPoly, new List<Corner>()
        {
            cornerC,
            cornerCA,
            cornerBC
        });
            newFaces[3] = new Face(tessPoly, new List<Corner>()
        {
            cornerAB,
            cornerBC,
            cornerCA
        });

            tessPoly.faces.AddRange(newFaces);

            // Yield if the time slice has been exceeded
            await yieldTimer.YieldOnTimeSlice();
            if (cancelRef.doCancel)
                return null;
        }

        await tessPoly.RecomputeEdgesAsync(cancelRef);
        return tessPoly;
        async UniTask<Corner> FindOrCreateCornerAsync(Polyhedron polyhedron, Vector3 vertex)
        {
          //  await UniTask.SwitchToThreadPool(); // Run this on a thread pool
            Corner corner = polyhedron.FindCornerByVertex(vertex);
            if (corner == null)
            {
                corner = new Corner(polyhedron, vertex);
                polyhedron.corners.Add(corner);
            }
            return corner;
        }
    }
    
    
    /// <summary>
    /// Sets all corner coordinates to have a length of 1 (from origin), without changing the direction from the origin.
    /// </summary>
    /// <param name="R"></param>
    public async UniTask SpherizeAsync(float R, CancelBoolRef cancelRef)
    {
        YieldTimer yieldTimer = new YieldTimer(cancelRef);
        foreach (Corner c in corners)
        {
            c.vertex = c.vertex.normalized * R;
            await yieldTimer.YieldOnTimeSlice();
        }
    }
    public void Spherize(float R)
    {
        foreach (Corner c in corners)
            c.vertex = c.vertex.normalized * R;
        /*
        // Optional smoothing step to reduce irregularities after spherizing
        int iterations = 5; // Number of smoothing passes
        for (int i = 0; i < iterations; i++)
        {
            foreach (Corner c in corners)
            {
                // Calculate the average position of neighboring vertices
                Vector3 avgPos = Vector3.zero;
                int neighborCount = c.touchingEdges.Count;

                for (int j = 0; j < neighborCount; j++)
                {
                    Vector3 neighborPos = c.touchingEdges[j].VectorFrom(c) + c.vertex;
                    avgPos += neighborPos;
                }

                if (neighborCount > 0)
                {
                    avgPos /= neighborCount;
                }

                // Blend towards the average to smooth out distortions
                c.vertex = Vector3.Lerp(c.vertex, avgPos.normalized * R, 0.1f);
            }
        }*/
    }

    //dual converts corners into faces, and faces into corners
    //Face to corner: the new corner will be "over" the center of the original face- but at the same distance from the polyhedron's center as original points
    // corner to face:  the face will have as many corners, as the original corner had edges touching it.
    //         the new face's corners will be those generated from the faces that used the original corner.
    public async UniTask<Polyhedron> DualAsync(CancelBoolRef cancelRef)
    {
        YieldTimer yieldTimer = new YieldTimer(cancelRef);
        Polyhedron dualPoly = new Polyhedron();
        List<Corner> newCorners = dualPoly.corners;//  this list- after first loop- we be the same size as- and in the same order as ``this.faces``
        List<Face> newFaces = dualPoly.faces;
        foreach (Face oldFace in this.faces)
        {

            List<Vector3> oldFaceCornerVerts = Corner.GetVerticies(oldFace.corners);
            Vector3 center = oldFaceCornerVerts.AvgPos().normalized;
            center *= oldFace.corners[0].vertex.magnitude;
            newCorners.Add(new Corner(dualPoly, center));
            // Yield if the time slice has been exceeded
            await yieldTimer.YieldOnTimeSlice();
        }

        foreach (Corner oldCorner in this.corners)
        {
            List<Corner> newFaceCorners = new List<Corner>();
            List<Face> oldCornersOldFaces = oldCorner.touchingFaces;
            //  Debug.Log("old corner at index: " + corners.IndexOf(oldCorner) + " touching " + oldCornersOldFaces.Count + "faces");

            foreach (Face oldFaceTouchingOldCorner in oldCornersOldFaces)
            {
                newFaceCorners.Add(newCorners[this.faces.IndexOf(oldFaceTouchingOldCorner)]);
            }
            newFaces.Add(new Face(dualPoly, newFaceCorners));
            // Yield if the time slice has been exceeded
            await yieldTimer.YieldOnTimeSlice();
        }
        foreach (Face f in newFaces)
        {
            f.ReOrderCornersClockWiseAroundCenterAndNormal();
            // Yield if the time slice has been exceeded
            await yieldTimer.YieldOnTimeSlice();
        }

        dualPoly.RecomputeEdges();
        //  dualPoly.CheckIntegrity();
        return dualPoly;
    }
    public Polyhedron Dual()
    {
        Polyhedron dualPoly = new Polyhedron();
        List<Corner> newCorners = dualPoly.corners;//  this list- after first loop- we be the same size as- and in the same order as ``this.faces``
        List<Face> newFaces = dualPoly.faces;
        foreach (Face oldFace in this.faces)
        {
           
            List<Vector3> oldFaceCornerVerts = Corner.GetVerticies(oldFace.corners);
            Vector3 center  = oldFaceCornerVerts.AvgPos().normalized;
            center *= oldFace.corners[0].vertex.magnitude;
            newCorners.Add(new Corner(dualPoly, center));
        }
       
        foreach (Corner oldCorner in this.corners)
        {
            List<Corner> newFaceCorners = new List<Corner>();
            List<Face> oldCornersOldFaces = oldCorner.touchingFaces;
          //  Debug.Log("old corner at index: " + corners.IndexOf(oldCorner) + " touching " + oldCornersOldFaces.Count + "faces");

            foreach (Face oldFaceTouchingOldCorner in oldCornersOldFaces)
            {
                newFaceCorners.Add(newCorners[this.faces.IndexOf(oldFaceTouchingOldCorner)]);
            }
            newFaces.Add(new Face(dualPoly, newFaceCorners));
        }
        foreach(Face f in newFaces)
            f.ReOrderCornersClockWiseAroundCenterAndNormal();

        dualPoly.RecomputeEdges();
      //  dualPoly.CheckIntegrity();
        return dualPoly;
    }
    public async UniTask<FacesAndNeighbors> GetFacesAndNeighborsAsync(CancelBoolRef cancelRef)
    {
        FacesAndNeighbors result = FacesAndNeighbors.CreateInstance<FacesAndNeighbors>();
        await ToMeshAsync(result, cancelRef);//ignore mesh output
        return result;
    }
    public FacesAndNeighbors GetFacesAndNeighbors()
    {
        FacesAndNeighbors result= FacesAndNeighbors.CreateInstance<FacesAndNeighbors>();
        ToMesh(result);//ignore mesh output
        return result;
    }

    /// <summary>
    /// Generate a unity mesh from this polyhedron.
    /// </summary>
    /// <param name="facesAndNeighbors">If provided, this reference will be populated with this polygon's face and neighbor info</param>
    /// <returns></returns>
    public async UniTask<Mesh> ToMeshAsync(FacesAndNeighbors facesAndNeighbors = null, CancelBoolRef cancelRef=null)
    {
        MeshData meshDataToUse = await ToMeshDataAsync(facesAndNeighbors, cancelRef);
        await UniTask.SwitchToMainThread();
        return meshDataToUse.ToMesh();
    }
    public Mesh ToMesh(FacesAndNeighbors facesAndNeighbors = null)
    {
        MeshData meshDataToUse = ToMeshData(facesAndNeighbors);
        return meshDataToUse.ToMesh();
    }
    /// <summary>
    /// Generate a unity mesh from this polyhedron.
    /// </summary>
    /// <param name="facesAndNeighbors">If provided, this reference will be populated with this polygon's face and neighbor info</param>
    /// <returns></returns>
    public async UniTask<MeshData> ToMeshDataAsync(FacesAndNeighbors facesAndNeighbors = null, CancelBoolRef cancelRef=null)
    {
        MeshData meshToUseRef = new MeshData();
        YieldTimer timer = new YieldTimer(cancelRef);
        bool doFacesAndNeighbors = (facesAndNeighbors != null);
        if (doFacesAndNeighbors)
            facesAndNeighbors.faceDetails = new List<FaceDetails>();
        List<Vector3> meshVerts = new List<Vector3>();
        List<int> meshTris = new List<int>();
        List<Vector3> meshNormals = new List<Vector3>();
        //     List<Vector2> meshUV0s = new List<Vector2>();
        List<Vector2> meshUV1s = new List<Vector2>();
        int faceCount = 0;
        //by face.... more than 3 corners we compute a "middle" point and connect multiple triangles to that
        foreach (Face currentFace in faces)
        {
            FaceDetails faceIndex = null;
            if (doFacesAndNeighbors)
            {
                faceIndex = new FaceDetails();
                faceIndex.index = faceCount++; //faces.IndexOf(currentFace);
                faceIndex.neighborIndices = new List<int>();
                foreach (Face neighbor in currentFace.neighbors)
                    faceIndex.neighborIndices.Add(faces.IndexOf(neighbor));
                faceIndex.normal = currentFace.normal;
                faceIndex.triangles = new List<int>(); // filled lower down
                facesAndNeighbors.faceDetails.Add(faceIndex);
            }

            //  currentFace.ReOrderCornersClockWiseAroundCenterAndNormal();
            int triIndexStart = meshVerts.Count;
            List<Vector3> faceVerts = Corner.GetVerticies(currentFace.corners);
            List<int> faceTris = new List<int>();
            Vector3 faceNormal = currentFace.normal;
            List<Vector3> faceNormals = new List<Vector3>();
            //List<Vector2> faceUV0s = new List<Vector2>();
            List<Vector2> faceUV1s = new List<Vector2>();


            if (currentFace.corners.Count != 3)
            {
                //compute center point of face (avg all points)
                Vector3 centerPt = faceVerts.AvgPos();
                // int centerPtIndex = faceVerts.Count;
                int numOriginalVerts = faceVerts.Count;
                //  faceVerts.Add(centerPt);
                List<Vector3> originalfaceVerts = faceVerts;
                faceVerts = new List<Vector3>();
                for (int i = 0; i < numOriginalVerts; i++)
                {
                    faceTris.Add(faceVerts.Count + triIndexStart);
                    faceVerts.Add(originalfaceVerts[(i + 0).CircularIndex(numOriginalVerts)]);

                    faceTris.Add(faceVerts.Count + triIndexStart);
                    faceVerts.Add(originalfaceVerts[(i + 1).CircularIndex(numOriginalVerts)]);
                    //faceTris.Add((i + 0).CircularIndex(numOriginalVerts) + triIndexStart);
                    //faceTris.Add((i + 1).CircularIndex(numOriginalVerts) + triIndexStart);

                    faceTris.Add(faceVerts.Count + triIndexStart);
                    faceVerts.Add(centerPt);
                }
            }
            else//3 pts- we can just add tri directly- no center needed
            {
                faceTris.Add(0 + triIndexStart);
                faceTris.Add(1 + triIndexStart);
                faceTris.Add(2 + triIndexStart);
            }

            if (doFacesAndNeighbors)
            {
                int triStartIndex = meshTris.Count;
                for (int i = 0; i < faceTris.Count; i += 3)
                    faceIndex.triangles.Add(i + triStartIndex);
            }

            Vector2 middleUV = Vector2.one * 0.5f;
            Vector3 faceCenter = currentFace.center;

            for (int i = 0; i < faceVerts.Count; i++)
            {
                faceNormals.Add(faceNormal);
                // faceUV0s.Add(Vector2.zero);// set later- by tri // faceVerts[i].SphericalUV());
                faceUV1s.Add((faceVerts[i].ProjectPointOntoPlane(faceNormal, faceCenter).normalized * 0.5f) + middleUV);
            }

            meshVerts.AddRange(faceVerts);
            meshTris.AddRange(faceTris);
            meshNormals.AddRange(faceNormals);
            //   meshUV0s.AddRange(faceUV0s);
            meshUV1s.AddRange(faceUV1s);
            await timer.YieldOnTimeSlice();
        }// end for each face


        meshToUseRef.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        meshToUseRef.meshVerts = meshVerts.ToArray();
        meshToUseRef.meshTris = meshTris.ToArray();
        meshToUseRef.meshNormals = meshNormals.ToArray();
        // newMesh.SetUVs(0, meshUV0s.ToArray());
        meshToUseRef.meshUV1s = meshUV1s.ToArray();
        meshToUseRef.meshUV2s = meshUV1s.ToArray();
        SphereizeMeshUV0s(meshToUseRef);
        meshToUseRef.name = faces.Count.ToString() + "FacePoly";
        if (doFacesAndNeighbors)
            meshToUseRef.facesAndNeighborsRef = facesAndNeighbors;

        


        //  TestCheckUv0Xs(newMesh);
        return meshToUseRef;
    }
    public MeshData ToMeshData(FacesAndNeighbors facesAndNeighbors = null)
    {
        MeshData meshToUseRef = new MeshData();
        bool doFacesAndNeighbors = (facesAndNeighbors != null);
        if (doFacesAndNeighbors)
            facesAndNeighbors.faceDetails = new List<FaceDetails>();
        List<Vector3> meshVerts = new List<Vector3>();
        List<int> meshTris = new List<int>();
        List<Vector3> meshNormals = new List<Vector3>();
   //     List<Vector2> meshUV0s = new List<Vector2>();
        List<Vector2> meshUV1s = new List<Vector2>();
        int faceCount = 0;
        //by face.... more than 3 corners we compute a "middle" point and connect multiple triangles to that
        foreach (Face currentFace in faces)
        {
            FaceDetails faceIndex=null;
            if (doFacesAndNeighbors)
            {
                faceIndex = new FaceDetails();
                faceIndex.index = faceCount++; //faces.IndexOf(currentFace);
                faceIndex.neighborIndices = new List<int>();
                foreach (Face neighbor in currentFace.neighbors)
                    faceIndex.neighborIndices.Add(faces.IndexOf(neighbor));
                faceIndex.normal = currentFace.normal;
                faceIndex.triangles = new List<int>(); // filled lower down
                facesAndNeighbors.faceDetails.Add(faceIndex);
            }

          //  currentFace.ReOrderCornersClockWiseAroundCenterAndNormal();
            int triIndexStart = meshVerts.Count;
            List<Vector3> faceVerts = Corner.GetVerticies(currentFace.corners);
            List<int> faceTris = new List<int>();
            Vector3 faceNormal = currentFace.normal;
            List<Vector3> faceNormals = new List<Vector3>();
            //List<Vector2> faceUV0s = new List<Vector2>();
            List<Vector2> faceUV1s = new List<Vector2>();


            if (currentFace.corners.Count != 3)
            {
                //compute center point of face (avg all points)
                Vector3 centerPt = faceVerts.AvgPos();
               // int centerPtIndex = faceVerts.Count;
                int numOriginalVerts = faceVerts.Count;
                //  faceVerts.Add(centerPt);
                List<Vector3> originalfaceVerts = faceVerts;
                faceVerts = new List<Vector3>();
                for (int i = 0; i < numOriginalVerts; i++)
                {
                    faceTris.Add(faceVerts.Count + triIndexStart);
                    faceVerts.Add(originalfaceVerts[(i + 0).CircularIndex(numOriginalVerts)]);
                    
                    faceTris.Add(faceVerts.Count + triIndexStart);
                    faceVerts.Add(originalfaceVerts[(i + 1).CircularIndex(numOriginalVerts)]);
                    //faceTris.Add((i + 0).CircularIndex(numOriginalVerts) + triIndexStart);
                    //faceTris.Add((i + 1).CircularIndex(numOriginalVerts) + triIndexStart);
                    
                    faceTris.Add(faceVerts.Count + triIndexStart);
                    faceVerts.Add(centerPt);
                }
            }
            else//3 pts- we can just add tri directly- no center needed
            { 
                faceTris.Add(0 + triIndexStart);
                faceTris.Add(1 + triIndexStart);
                faceTris.Add(2 + triIndexStart);
            }

            if (doFacesAndNeighbors)
            {
                int triStartIndex = meshTris.Count;
                for (int i = 0; i < faceTris.Count; i += 3)
                    faceIndex.triangles.Add(i + triStartIndex);
            }

            Vector2 middleUV = Vector2.one * 0.5f;
            Vector3 faceCenter = currentFace.center;
            
            for (int i = 0; i < faceVerts.Count; i++)
            {
                faceNormals.Add(faceNormal);
               // faceUV0s.Add(Vector2.zero);// set later- by tri // faceVerts[i].SphericalUV());
                faceUV1s.Add((faceVerts[i].ProjectPointOntoPlane(faceNormal, faceCenter).normalized * 0.5f) + middleUV);
            }

            meshVerts.AddRange(faceVerts);
            meshTris.AddRange(faceTris);
            meshNormals.AddRange(faceNormals);
         //   meshUV0s.AddRange(faceUV0s);
            meshUV1s.AddRange(faceUV1s);
        }// end for each face


        meshToUseRef.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        meshToUseRef.meshVerts = meshVerts.ToArray();
        meshToUseRef.meshTris = meshTris.ToArray();
        meshToUseRef.meshNormals = meshNormals.ToArray();
        // newMesh.SetUVs(0, meshUV0s.ToArray());
        meshToUseRef.meshUV1s = meshUV1s.ToArray();
        meshToUseRef.meshUV2s = meshUV1s.ToArray();
        SphereizeMeshUV0s(meshToUseRef);
        meshToUseRef.name = faces.Count.ToString() + "FacePoly";
        if (doFacesAndNeighbors)
            meshToUseRef.facesAndNeighborsRef = facesAndNeighbors;
            


      //  TestCheckUv0Xs(newMesh);
        return meshToUseRef;
    }
    

    static public void SphereizeMeshUV0sWithTests(Mesh inputMesh)
    {
        int[] triangleArray = inputMesh.triangles;
        Vector3[] vertexArray = inputMesh.vertices;
        Vector2[] uvArray = new Vector2[vertexArray.Length];

        for (int i = 0; i < triangleArray.Length; i += 3)
        {
            int index0 = triangleArray[i];
            int index1 = triangleArray[i + 1];
            int index2 = triangleArray[i + 2];

            //spherical uv's wrap X from 0 to 1 around the planet horizontally (around the equator)
            // spherical uv'sn have Y =1 at north pole, y=0 at south pole, and y=0.5 at equator             
            Vector2 uv0 = vertexArray[index0].CylindricalUV();
            Vector2 uv1 = vertexArray[index1].CylindricalUV();
            Vector2 uv2 = vertexArray[index2].CylindricalUV();

            //account for triangles that cross the 0,1 meridian : assumes texture sample mode is "repeat"
            float UV0xNew = uv0.x;
            float UV1xNew = uv1.x;
            float UV2xNew = uv2.x;

            if (uv1.x - uv0.x > 0.5f)
                UV1xNew -= 1f;
            if (uv1.x - uv0.x < -0.5f)
                UV1xNew += 1f;

            if (uv2.x - uv0.x > 0.5f)
                UV2xNew -= 1f;
            if (uv2.x - uv0.x < -0.5f)
                UV2xNew += 1f;
            /*
            if (uv2.x - uv1.x > 0.5f)
                UV2xNew -= 1f;
            if (uv2.x - uv1.x < -0.5f)
                UV2xNew += 1f;
            */

            #region troubleShooting
            bool check = false;
            if (Mathf.Abs(uv1.x - uv0.x) > 0.5f)
                check = true;
            if (Mathf.Abs(uv2.x - uv0.x) > 0.5f)
                check = true;
            if (Mathf.Abs(uv2.x - uv1.x) > 0.5f)
                check = true;

           
            bool stillOff = false;
            if (Mathf.Abs(UV1xNew - UV0xNew) > 0.5f)
            {
                stillOff = true;
                check = true;
            }
            if (Mathf.Abs(UV1xNew - UV2xNew) > 0.5f)
            {
                stillOff = true;
                check = true;
            }
            if (Mathf.Abs(UV2xNew - UV0xNew) > 0.5f)
            {
                stillOff = true;
                check = true;
            }

            if(check)
            {
                string s = "Modified uvs for tri index: " + (i/3);
                s+="\nBEFORE Tri UVx's: " + uv0.x.ToString("0.00")
                    + "," + uv1.x.ToString("0.00") + "," + uv2.x.ToString("0.00");
                s+="\nAfter Tri UVx's: " + UV0xNew.ToString("0.00")
                    + "," + UV1xNew.ToString("0.00") + "," + UV2xNew.ToString("0.00");
                s += "\nStill spans 1/2 rotation: " + stillOff;
                s += "\nVertIndexes: " + index0 +","+ index1 + "," + index2;
                    
                Debug.Log(s);
            }

            #endregion
            uv0.x = UV0xNew;
            uv1.x = UV1xNew;
            uv2.x = UV2xNew;
            uv0 = new Vector2(UV0xNew, uv0.y);
            uv1 = new Vector2(UV1xNew, uv1.y);
            uv2 = new Vector2(UV2xNew, uv2.y);
            
            if (uv0.y.CloseEqual(1) || uv0.y.CloseEqual(0))
                uv0.x = (uv1.x + uv2.x) * 0.5f; //avg
            if (uv1.y.CloseEqual(1) || uv1.y.CloseEqual(0))
                uv1.x = (uv0.x + uv2.x) * 0.5f; //avg
            if (uv2.y.CloseEqual(1) || uv2.y.CloseEqual(0))
                uv2.x = (uv1.x + uv0.x) * 0.5f; //avg
            

            uvArray[index0] = uv0;
            uvArray[index1] = uv1;
            uvArray[index2] = uv2;
        }
        inputMesh.SetUVs(0, uvArray,0,uvArray.Length );

        List<Vector2> uvCheckArray = new List<Vector2>();
        inputMesh.GetUVs(0, uvCheckArray);//.uv;
        for (int i = 0; i < uvCheckArray.Count; i++)
            if (uvCheckArray[i] != uvArray[i])
                Debug.Log("SanityFailure");
    }
    static public void SphereizeMeshUV0s(MeshData inputMesh)
    {
        int[] triangleArray = inputMesh.meshTris;
        Vector3[] vertexArray = inputMesh.meshVerts;
        Vector2[] uvArray = new Vector2[vertexArray.Length];

        for (int i = 0; i < triangleArray.Length; i += 3)
        {
            int index0 = triangleArray[i];
            int index1 = triangleArray[i + 1];
            int index2 = triangleArray[i + 2];

            //spherical uv's wrap X from 0 to 1 around the planet horizontally (around the equator)
            // spherical uv'sn have Y =1 at north pole, y=0 at south pole, and y=0.5 at equator             
            Vector2 uv0 = vertexArray[index0].CylindricalUV();
            Vector2 uv1 = vertexArray[index1].CylindricalUV();
            Vector2 uv2 = vertexArray[index2].CylindricalUV();



            bool closeToPole = false;
            if (uv0.y.CloseEqual(1) || uv0.y.CloseEqual(0))
            {
                uv0.x = (uv1.x + uv2.x) * 0.5f; //avg
                closeToPole = true;
            }
            if (uv1.y.CloseEqual(1) || uv1.y.CloseEqual(0))
            {
                uv1.x = (uv0.x + uv2.x) * 0.5f; //avg
                closeToPole = true;
            }
            if (uv2.y.CloseEqual(1) || uv2.y.CloseEqual(0))
            {
                uv2.x = (uv1.x + uv0.x) * 0.5f; //avg
                closeToPole = true;
            }

            if (!closeToPole)
            {
                //account for triangles that cross the 0,1 meridian : assumes texture sample mode is "repeat"

                if (uv1.x - uv0.x > 0.5f)
                    uv1.x -= 1f;
                else if (uv1.x - uv0.x < -0.5f)
                    uv1.x += 1f;

                if (uv2.x - uv0.x > 0.5f)
                    uv2.x -= 1f;
                else if (uv2.x - uv0.x < -0.5f)
                    uv2.x += 1f;
            }

            uvArray[index0] = uv0;
            uvArray[index1] = uv1;
            uvArray[index2] = uv2;
        }
        inputMesh.meshUV0s= uvArray;

    }

    static void TestCheckUv0Xs(Mesh mesh)
    {
        int[] triangleArray = mesh.triangles;
        List<Vector2> uvArray = new List<Vector2>();
        mesh.GetUVs(0, uvArray);//.uv;
        for (int i = 0; i < triangleArray.Length; i += 3)
        {
            int index0 = triangleArray[i];
            int index1 = triangleArray[i + 1];
            int index2 = triangleArray[i + 2];
            Vector2 uv0 = uvArray[index0];
            Vector2 uv1 = uvArray[index1];
            Vector2 uv2 = uvArray[index2];
            bool failcheck = false;
            if (Mathf.Abs(uv1.x - uv0.x) > 0.5f)
                failcheck = true;
            if (Mathf.Abs(uv2.x - uv0.x) > 0.5f)
                failcheck = true;
            if (Mathf.Abs(uv2.x - uv1.x) > 0.5f)
                failcheck = true;
            if(failcheck)
                Debug.Log("FAILCHECK Tri " +(i/3) +"   UVx's: " + uv0.x.ToString("0.00")
                    + "," + uv1.x.ToString("0.00") + "," + uv2.x.ToString("0.00") +
                    "\nVertIndexes: " + index0 +", "+ index1 + ", " + index2);
        }


    }

    const float radToRots = 1.0f / (Mathf.PI * 2.0f);

    public bool CheckIntegrity(Transform debugDragAt, Gradient orderGradient)
    {
        Debug.Log("Checking Integrity");
        if (corners.Count < 3)
        {
            Debug.Log("Poly has less than 3 corners- invalid");
            return false;
        }
        if (corners.Count == 3)
        {
            Debug.Log("Poly has exactly 3 corners- valid for 2d objects only");
        }
        foreach (Corner c in corners)
        {
            if (c.touchingFaces.Count < 1)
            {
                Debug.Log("Found corner touching no faces: invalid");
                return false;
            }
            if (c.touchingFaces.Count < 3)
            {
                Debug.Log("Found corner touching less than 3 faces- valid for 2d objects only");
            }

            if (c.touchingEdges.Count < 2)
            {
                Debug.Log("Found corner touching one or no edges: invalid");
                return false;
            }
            else
            {
                List<Edge> edgesTouchingCorner = c.touchingEdges;
                if (edgesTouchingCorner.Count < 3)
                {
                    Debug.Log("Found corner touching less than 3 edges- valid for 2d objects only");
                }

                Debug.Log("CheckingEdges around corner");
                //float angle = Vector3.SignedAngle(edgesTouchingCorner[0].VectorFrom(c), edgesTouchingCorner[1].VectorFrom(c),c.vertex);
                float angle = Vector3.Angle(edgesTouchingCorner[0].VectorFrom(c), edgesTouchingCorner[1].VectorFrom(c));
                //Debug.Log("    AngleBetween edge 0: " + edgesTouchingCorner[0].VectorFrom(c) + " and edge 1: " + edgesTouchingCorner[1].VectorFrom(c) + " is " + angle);
                for (int i = 0; i < edgesTouchingCorner.Count; i++)
                {
                    int j = i + 1;
                    if(j >= edgesTouchingCorner.Count) j = 0;
                    
                    //float this2EdgeAngle = Vector3.SignedAngle(edgesTouchingCorner[i].VectorFrom(c), edgesTouchingCorner[j].VectorFrom(c), c.vertex);
                    float this2EdgeAngle = Vector3.Angle(edgesTouchingCorner[i].VectorFrom(c), edgesTouchingCorner[j].VectorFrom(c));
                    //Debug.Log("    AngleBetween edge " + i + ": "+ edgesTouchingCorner[i].VectorFrom(c) + " and edge " + j + ": "+ edgesTouchingCorner[j].VectorFrom(c) + " is " + this2EdgeAngle);
                    float angleDiff = Mathf.Abs(angle - this2EdgeAngle);
                    if (angleDiff > 0.001f)//angle != Vector3.Angle(c.touchingEdges[i].VectorFrom(c), c.touchingEdges[j].VectorFrom(c)))
                    {
                        Debug.Log("    Found Corner with touching edges that are not evenly angled around the corner ("+angle+","+ this2EdgeAngle + ")");
                        int count = 0;
                        for (int x = 0; x < edgesTouchingCorner.Count; x++)
                        //    foreach (Edge eachEdge in edgesTouchingCorner)
                        {
                            Edge eachEdge = edgesTouchingCorner[x];
                            float frac = ((float)count++ / (float)edgesTouchingCorner.Count);
                            Color col = orderGradient.Evaluate(frac);//   Color.green *  (0.5f+(count++/ edgesTouchingCorner.Count*2));
                            Debug.DrawLine(debugDragAt.TransformPoint(eachEdge.endpointA.vertex),
                                    debugDragAt.TransformPoint(eachEdge.endpointB.vertex), col, 500);
                            
                            int y = x + 1;
                            if (y >= edgesTouchingCorner.Count) y = 0;
                            //float aAngle = Vector3.SignedAngle(edgesTouchingCorner[x].VectorFrom(c), edgesTouchingCorner[y].VectorFrom(c), c.vertex);
                            float aAngle = Vector3.Angle(edgesTouchingCorner[x].VectorFrom(c), edgesTouchingCorner[y].VectorFrom(c));
                          //  Debug.Log("drawing edge :" + eachEdge.endpointA.vertex + " to " + eachEdge.endpointB.vertex + " Angle to next: "+aAngle);
                        }
                        Debug.Log("just drew " + count + " edges");
                       // return false;
                        break;
                    }
                }

            }
        }

        int faceIndex = 0;
        foreach (Face f in faces)
        {
            List<Edge> faceEdges = f.edges;
            if (faceEdges.Count < 3)
                Debug.Log("Found face with less than 3 edges- invalid.");
            float len = faceEdges[0].length;
            for (int i = 1; i < faceEdges.Count; i++)
            {
                float checkLen = faceEdges[i].length;
                if (Mathf.Abs(checkLen- len) > .0001f* len)
                {
                  //  Debug.Log("Found Face at index " + faceIndex + ", with different length Edges("+len+","+ checkLen + "): Non-Uniform");
                    foreach (Edge eachEdge in faceEdges)
                    {
                        Color col = Color.white;
                        if (eachEdge == faceEdges[0])
                            col = Color.green;
                        if (eachEdge.length != len)
                            col = Color.red;
                        Debug.DrawLine(debugDragAt.TransformPoint(eachEdge.endpointA.vertex),
                                debugDragAt.TransformPoint(eachEdge.endpointB.vertex), col, 500);

                    }

                    break;
                }
            }


            faceIndex++;
        }

        float lenSum = 0;
        int ecount = 0;
        float maxLen = 0;
        float minLen = float.MaxValue;
        foreach (Edge e in edges)
        {
            if (e.touchingFaces.Count < 2)
                Debug.Log("Found edge touching less than 2 faces- invalid for 3d objects");
            float len = e.length;
            lenSum += len;
            ecount++;
            if (len > maxLen) maxLen = len;
            if (len < minLen) minLen = len;
        }
        float avglen = lenSum / (float)ecount;
        float maxDiff = Mathf.Abs(maxLen - minLen);
        Debug.Log("All poly Edges-  avg len: " + avglen + " min len: " + minLen + "  max len: " + maxLen + "  uniform: "+ (maxDiff<.0001f));
        return true;
    }

    public static Polyhedron ComputeIcosahedron()
    {
        //  φ = (1+√5)/2
        float φ = (1 + (float)Mathf.Sqrt(5)) / 2;
        float neg_a = φ * -1.0f;
        //Ico icos = new Ico();
        //(0, ±1, ±φ)
        //(±1, ±φ, 0)
        //(±φ, 0, ±1)
        // ModelMesh ico_mesh=ico.Meshes[0];

        Vector3[] vertices = new Vector3[12];
        Vector3[] normals = new Vector3[12];
        Vector2[] uv = new Vector2[12];

        //Initialize the custom vertex values for the triangle strip.
        //float[5] possible_values={neg_a,-1,0,1,a};
        //            Vector3 vert=new Vector3(0.0f, -1.0f, neg_a)
        //          triangleStripVertices[0] = new VertexPositionNormalColored(vert,Color.Red,vert);

        //icoshedron has 12 VertexBuffer, and 20 faces

        int vert_counter = 0;

        //golden rect on Y axis XZ PLANE
        for (float z = neg_a; z <= φ; z += 2.0f * φ)
        {
            for (float x = -1.0f; x <= 1.0f; x += 2.0f)
            {
                Vector3 vert = new Vector3(x, 0.0f, z);
                Vector3 vert_normal = new Vector3(x, 0.0f, z);
                vertices[vert_counter++] = vert;// new VertexPositionNormalColoredTextured(vert, Color.Red, vert_normal);
            }
        }

        //golden rect on Z axis YX PLANE
        for (float x = neg_a; x <= φ; x += 2.0f * φ)
        {
            for (float y = -1.0f; y <= 1.0f; y += 2.0f)
            {
                Vector3 vert = new Vector3(x, y, 0.0f);
                Vector3 vert_normal = new Vector3(x, y, 0.0f);
                vertices[vert_counter++] = vert;// new VertexPositionNormalColoredTextured(vert, Color.Red, vert_normal);
            }
        }
        //golden rect on X axis YZ PLANE
        for (float y = neg_a; y <= φ; y += 2.0f * φ)
        {
            for (float z = -1.0f; z <= 1.0f; z += 2.0f)
            {
                Vector3 vert = new Vector3(0.0f, y, z);
                Vector3 vert_normal = new Vector3(0.0f, y, z);
                vertices[vert_counter++] = vert;// new VertexPositionNormalColoredTextured(vert, Color.Red, vert_normal);
            }
        }

        //normalize al the "normals" (set magnitude = to 1- which means its a "dirction only vector")
        for (int i = 0; i < vertices.Length; i++)
        {
            vertices[i].Normalize();
            normals[i] = vertices[i];
            float u = 0.5f + (Mathf.Atan2(vertices[i].x, vertices[i].z) / (2 * Mathf.PI));
            float v = 0.5f - (Mathf.Asin(vertices[i].y) / (Mathf.PI));
            uv[i] = new Vector2(u, v);
        }
        int[] triangles = new int[20 * 3] {
                0, 4, 8,
                4, 9, 8,
                8, 9, 6,
                9, 3, 6,
                6, 3, 7,
                1, 6, 7,
                8, 6, 1,
                0, 8, 1,
                0, 1, 10,
                1, 7, 10,
                10, 7, 11,
                7, 3, 11,
                11, 3, 2,
                5, 11, 2,
                10, 11, 5,
                0, 10, 5,
                0, 5, 4,
                5, 2, 4,
                4, 2, 9,
                2, 3, 9 };

        //oops, above stripindex list is drawing the inside of the shape- this loop transposes stripindecies to the outside
        for (int i = 0; i < 20 * 3; i = i + 3)
        {
            int a = triangles[i];
            triangles[i] = triangles[i + 2];
            triangles[i + 2] = a;
        }


        //for poly- convert each tri to it's own list
        List<List<int>> facesByIndex = new List<List<int>>();
        for (int i = 0; i < 20; i++)
        {
            List<int> triList = new List<int>()
            {
                triangles[i*3+0],
                triangles[i*3+1],
                triangles[i*3+2]
            };
            facesByIndex.Add(triList);
        }
        return new Polyhedron(new List<Vector3>(vertices), facesByIndex);

    }
    public static Polyhedron ComputeFlatSquare()
    {
        return new Polyhedron(new List<Vector3>()
                { Vector3.zero, Vector3.right, new Vector3(1,1,0), Vector3.up },
                new List<List<int>>() { new List<int> { 0, 2, 1, 3 } }); //out of order test
    }
    public void OnAfterDeserialize()
    {
        //circular references bad for serialization do it now
        foreach (Face f in faces)
            f.SetPoly(this);
        foreach (Edge f in edges)
            f.SetPoly(this);
        foreach (Corner f in corners)
            f.SetPoly(this);
    }
    public void OnBeforeSerialize() { }
}
