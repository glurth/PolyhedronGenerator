using System.Collections;
//using System.Collections.Generic;
using UnityEngine;
using System.Threading;
//using System.Threading.Tasks;
using System;
using Cysharp.Threading.Tasks;
using EyE.UnityAssetTypes;

public class PolyhedronMaker : MonoBehaviour
{
    public Mesh startingMesh;//optional: will be used during start to compute a starting polyhedron, IF startingPoly is null.
    public PolyhedronSO startingPoly;//optional: will be used as the starting polyhedron, if provided.

    public bool defaultToTetraNotIcosa;// if neither of the above two are selected- adefault platonic shape will be usedas the starting polyhedron
    
    public Polyhedron poly;
    public MeshFilter mf;
    public FacesAndNeighbors facesAndNeighborsOnMesh;
    public int numberOfIteration = 3;
    [Range(0,0.5f)]
    public float fraction = 1f / 3f;


    bool isProcessing = false;
    async UniTask ProcessStart()
    {
        await UniTask.SwitchToMainThread();
        isProcessing = true;
        await UniTask.SwitchToTaskPool();
    }
    async UniTask ProcessDone()
    {
      //  processedMeshData = await poly.ToMeshDataAsync(facesAndNeighborsOnMesh, cancelRef);
        await UniTask.SwitchToMainThread();
        isProcessing = false;
    }


    // Start is called before the first frame update
    void Start()
    {

        facesAndNeighborsOnMesh = FacesAndNeighbors.CreateInstance<FacesAndNeighbors>();
        if (startingPoly != null)
        {
            poly = new Polyhedron(startingPoly.poly);
            SetDrawnMesh(startingPoly.mesh);
            facesAndNeighborsOnMesh = startingPoly.facesAndNeighbors;
        }
        else
        {
            if (startingMesh != null)
            {
                poly = new Polyhedron(startingMesh);// startingPoly;
                SetDrawnMesh(startingMesh);
            }
            else
            {
                if (defaultToTetraNotIcosa)
                    poly = Polyhedron.ComputeTetrahedron();
                else
                    poly = Polyhedron.ComputeIcosahedron();
                SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
            }
        }
       
        lastFraction = fraction;
    }
    /*
    public void GenerateTessAndDualFromStart()
    {
        StartCoroutine(DoTessAndDuals());
    }
    public void GenerateTruncAndDualFromStart()
    {
        StartCoroutine(DoTruncateAndDuals());
    }
    public void GenerateAlternatingTruncOrTessAndDualFromStart()
    {
        StartCoroutine(DoAlternatingTruncateAndTesselate(true));
    }

    public void GenerateAlternatingTessOrTruncAndDualFromStart()
    {
        StartCoroutine(DoAlternatingTruncateAndTesselate(false));
    }

    public IEnumerator DoTessAndDuals()
    {
       
        //poly = startingPoly;
        poly = poly.Dual();
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        SaveSOs();
        for (int i = 0; i < autoInterations; i++)
        {
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            poly = poly.TesselateTriangleByEdgeMiddles();
            yield return new WaitForSeconds(.1f);
            poly.Spherize(1f);
            yield return new WaitForSeconds(.1f);
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
            SaveSOs();
            yield return new WaitForSeconds(.1f);
        }

    }
    public IEnumerator DoTruncateAndDuals()
    {

       // poly = startingPoly;//icosa
        poly = poly.Dual();//dodeca
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        SaveSOs();

        for (int i = 0; i < autoInterations; i++)
        {
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            poly = poly.Trunc();
            yield return new WaitForSeconds(.1f);
            SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
            SaveSOs();
            yield return new WaitForSeconds(.1f);
        }

    }

    */
    private int counter = 0;
    private bool threadCompleted = false;
    
    // Generalized coroutine to run any function in a background thread with a timer
    private System.Collections.IEnumerator RunWithBackgroundTaskCoroutine(Action backgroundTask, int timerSeconds)
    {
        // Reset completion flag
        threadCompleted = false;

        // Start the specified function in a background thread
        Thread backgroundThread = new Thread(() =>
        {
            backgroundTask(); // Run the passed-in function
            threadCompleted = true; // Set flag to true when function completes
        });
        backgroundThread.Start();

        while(true)
        {
            // Increment the counter each second
            counter++;
            Debug.Log($"Processing Tick: {counter}");

            // Reduce the countdown timer
            yield return new WaitForSeconds(timerSeconds);

            // Exit if the thread has completed
            if (threadCompleted)
            {
                Debug.Log("Background task completed.");
                break;
            }
        }

        // Wait until the background thread completes if it didn't finish within the timer
        while (!threadCompleted)
        {
            yield return null; // Yield control back to Unity each frame
        }

        Debug.Log("Coroutine finished.");
    }

    /*public IEnumerator DoTriEdgeTessThenSpherizeOnceThenDualOnce()
    {

       // poly = startingPoly;
        for (int i = 0; i < numberOfIteration; i++)
        {
            yield return RunWithBackgroundTaskCoroutine(() =>{poly = poly.TesselateTriangleByEdgeMiddles();}, 3);
            yield return new WaitForSeconds(.1f);
        }
        yield return new WaitForSeconds(.1f);
        yield return RunWithBackgroundTaskCoroutine(() => { poly.Spherize(1f); }, 3);
        yield return new WaitForSeconds(.1f);
        yield return RunWithBackgroundTaskCoroutine(() => {poly = poly.Dual(); }, 3);
        yield return new WaitForSeconds(.1f);
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        //Mesh output=null;
        //yield return RunWithBackgroundTaskCoroutine(() => { output = poly.ToMesh(facesAndNeighborsOnMesh); }, 3);
        //SetDrawnMesh(output);
        SaveSOs();
    }*/
    public IEnumerator DoTessThenSpherizeOnce()
    {
        isProcessing = true;
        Debug.Log("Starting tesselation");
       // poly = startingPoly;
        for (int i = 0; i < numberOfIteration; i++)
        {
            Debug.Log("Starting tess pass ["+i+"]");
            yield return RunWithBackgroundTaskCoroutine(() => { poly = poly.TesselateTriangleByEdgeMiddles(); }, 3);
            yield return new WaitForSeconds(.1f);
        }
        yield return new WaitForSeconds(.1f);
        yield return RunWithBackgroundTaskCoroutine(() =>{poly.Spherize(1f);}, 3);
        yield return new WaitForSeconds(.1f);
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        SaveSOs();
        isProcessing = false;
    }
    public IEnumerator DoAlternatingTruncateAndTesselate(bool toggle = true)
    {
        isProcessing = true;
        //  poly = startingPoly;//icosa
        poly = poly.Dual();//dodeca
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        SaveSOs();
        
        for (int i = 0; i < numberOfIteration; i++)
        {
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            if (toggle)
                poly = poly.Trunc();
            else
            {
                poly = poly.TesselateTriangleByEdgeMiddles();
                yield return new WaitForSeconds(.1f);
                poly.Spherize(1f);
                yield return new WaitForSeconds(.1f);
                poly = poly.Dual();
            }
            toggle = !toggle;
            yield return new WaitForSeconds(.1f);
            SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
            SaveSOs();
            yield return new WaitForSeconds(.1f);
        }
        isProcessing = false;
    }
    public IEnumerator DoMultipleTriEdgeTessSpherizeDual()
    {
        isProcessing = true;
        poly = poly.Dual();
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        SaveSOs();
        for (int i = 0; i < numberOfIteration; i++)
        {
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            poly = poly.TesselateTriangleByEdgeMiddles();
            yield return new WaitForSeconds(.1f);
            poly.Spherize(1f);
            yield return new WaitForSeconds(.1f);
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
            SaveSOs();
            yield return new WaitForSeconds(.1f);
        }
        //processedMeshData = poly.ToMeshData(facesAndNeighborsOnMesh);
        isProcessing = false;
    }
    public IEnumerator DoRadialFaceTessSpherizeDual()
    {
        isProcessing = true;
        poly = poly.Dual();
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        SaveSOs();
        for (int i = 0; i < numberOfIteration; i++)
        {
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            poly = poly.TesselateFacesRadial();
            yield return new WaitForSeconds(.1f);
            poly.Spherize(1f);
            yield return new WaitForSeconds(.1f);
            poly = poly.Dual();
            yield return new WaitForSeconds(.1f);
            SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
            SaveSOs();
            yield return new WaitForSeconds(.1f);
        }
        isProcessing = false;
    }
    float lastFraction = -1f;

    // Update is called once per frame
    bool wasProcessingLastFrame = false;
    void Update()
    {

        if (wasProcessingLastFrame && !isProcessing && processedMeshData!=null)
        {
            SetDrawnMesh(processedMeshData.ToMesh());
        }
        wasProcessingLastFrame = isProcessing;
        
    }

    
    MeshData processedMeshData;

    Polyhedron.CancelBoolRef cancelRef = new Polyhedron.CancelBoolRef();
    // Method to cancel the task
    void CancelTask()
    {
        cancelRef.doCancel = true;
        cancelRef = new Polyhedron.CancelBoolRef();
    }

    public async UniTask TruncatePolyAsync()
    {
        poly = await poly.TruncAsync(fraction,cancelRef);// 1f / Mathf.Tan(Mathf.Deg2Rad * 360/5));
                                    //SetDrawnMesh(
        processedMeshData = poly.ToMeshData(facesAndNeighborsOnMesh);//);
    }
    public void TruncatePoly()
    {
        poly = poly.Trunc(fraction);// 1f / Mathf.Tan(Mathf.Deg2Rad * 360/5));
                                    //SetDrawnMesh(
        processedMeshData = poly.ToMeshData(facesAndNeighborsOnMesh);//);
    }


    public async UniTask DualPolyAsync()
    {
        await ProcessStart();
        poly = await poly.DualAsync(cancelRef);
        await ProcessDone();
    }
    public void DualPoly()
    {
        poly = poly.Dual();
        processedMeshData = poly.ToMeshData(facesAndNeighborsOnMesh);
   //     SetDrawnMesh(m);
    }


    public async UniTask TriangleTessalateAsync()
    {
        await ProcessStart();
        poly = await poly.TesselateTriangleByEdgeMiddlesAsync(cancelRef);
        await ProcessDone();
    }
    public void TriangleTessalate()
    {
        poly = poly.TesselateTriangleByEdgeMiddles();//.TesselateTriangles();
        //SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
        processedMeshData = poly.ToMeshData(facesAndNeighborsOnMesh);
    }
    
    
    public async void RadialFaceTessalateAsync()
    {
        await ProcessStart();
        poly = await poly.TesselateFacesRadialAsync(cancelRef);
        await ProcessDone();
    }
    public void RadialFaceTessalate()
    {
        poly = poly.TesselateFacesRadial();//.TesselateTriangles();
        poly.ToMeshData(facesAndNeighborsOnMesh);
        //        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
    }

    public async UniTask SpherizeAsync()
    {
        await ProcessStart();
        await poly.SpherizeAsync(1f,cancelRef);
        await ProcessDone();
    }

    public void Spherize()
    {
        poly.Spherize(1f);
        processedMeshData = poly.ToMeshData(facesAndNeighborsOnMesh);
        //SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
    }

    public async UniTask RecomputePolyMeshAsync()
    {
        await ProcessStart();
        SetDrawnMesh(await poly.ToMeshAsync(facesAndNeighborsOnMesh));
        await ProcessDone();
    }
    public void RecomputePolyMesh()
    {
        SetDrawnMesh(poly.ToMesh(facesAndNeighborsOnMesh));
    }

    public Gradient orderGradient = new Gradient();
    public void CheckPoly()
    {
        Debug.Log("Checking poly");
        poly.CheckIntegrity(transform, orderGradient);
    }

    void SetDrawnMesh(Mesh mesh)
    {
        mf.sharedMesh = mesh;
        MeshCollider mc;
        if (TryGetComponent<MeshCollider>(out mc))
            mc.sharedMesh = mesh;
    }

    public void SaveSOs()
    {
        string namePrefix = poly.faces.Count.ToString();
        if (mf.sharedMesh == null)
        {
            Debug.LogWarning("Unable to save mesh asset: unable to find one assigned to meshFilter.sharedMesh");
        }
        else
        {
            UnityEditor.AssetDatabase.CreateAsset(mf.sharedMesh, "assets/" + namePrefix + "Mesh.asset");
        }
       FacesAndNeighbors savedInstance = FacesAndNeighbors.Instantiate<FacesAndNeighbors>(facesAndNeighborsOnMesh);
       savedInstance.meshRef = mf.sharedMesh;
       UnityEditor.AssetDatabase.CreateAsset(savedInstance, "Assets/" + namePrefix + "FacesAndNeighbors.asset");

        PolyhedronSO polyhedronSO = PolyhedronSO.CreateInstance<PolyhedronSO>();
        polyhedronSO.poly = poly;
        polyhedronSO.facesAndNeighbors = savedInstance;
        polyhedronSO.mesh = mf.sharedMesh;
       UnityEditor.AssetDatabase.CreateAsset(polyhedronSO, "Assets/" + namePrefix + "Polyhedron.asset");

    }

    void OnGUI()
    {
        GUI.enabled = !isProcessing;
        /*
        // Generate button for each specific method, explicitly defined
        if (GUILayout.Button("GenerateTessAndDualFromStart"))
        {
            GenerateTessAndDualFromStart();
        }
        if (GUILayout.Button("GenerateTruncAndDualFromStart"))
        {
            GenerateTruncAndDualFromStart();
        }
        if (GUILayout.Button("GenerateAlternatingTruncOrTessAndDualFromStart"))
        {
            GenerateAlternatingTruncOrTessAndDualFromStart();
        }
        if (GUILayout.Button("GenerateAlternatingTessOrTruncAndDualFromStart"))
        {
            GenerateAlternatingTessOrTruncAndDualFromStart();
        }
        */
        GUILayout.Label($"Number of Iterations: {numberOfIteration}");
        numberOfIteration = Mathf.RoundToInt(GUILayout.HorizontalSlider(numberOfIteration, 1, 10));

        if (GUILayout.Button(new GUIContent("Multiple TriangleTesselation, Spherize, Dual Once","Use on polyhedrons with triangular faces")))
            StartCoroutine(DoMultipleTriEdgeTessSpherizeDual());
        if (GUILayout.Button(new GUIContent("Multiple RadialFaceTesselation, Spherize, Dual Once", "Use on polyhedrons with triangular faces")))
            StartCoroutine(DoRadialFaceTessSpherizeDual());
        

        if (GUILayout.Button(new GUIContent("Multiple TriangleTesselation then, Spherize Once", "Use on polyhedrons with triangular faces")))
            StartCoroutine(DoTessThenSpherizeOnce());
        
        if (GUILayout.Button("TruncatePoly"))
        {
            TruncatePoly();
            //LaunchAsync(TruncatePoly);
        }
        if (GUILayout.Button("DualPoly"))
        {
            DualPolyAsync();
            //LaunchAsync(DualPoly);
        }
        if (GUILayout.Button("Triangle Edge Tesselate"))
        {
            TriangleTessalateAsync();
           // LaunchAsync(TriangleTessalate);
        }
        if (GUILayout.Button("Radial Face Tesselate"))
        {
            RadialFaceTessalateAsync();
            //LaunchAsync(RadialFaceTessalate);
        }

        if (GUILayout.Button("Spherize"))
        {
            SpherizeAsync();
            //LaunchAsync(Spherize);
        }
        if (GUILayout.Button("Smooth"))
        {
            poly.Smooth();
            //LaunchAsync(Spherize);
        }
        if (GUILayout.Button("IcoShpere"))
        {
            Mesh ico=poly.Icosphere(numberOfIteration);
            SetDrawnMesh(ico);
            //LaunchAsync(Spherize);
        }

        if (GUILayout.Button("RecomputePolyMesh"))
        {
            RecomputePolyMesh();
        }
        if (GUILayout.Button("CheckPoly"))
        {
            CheckPoly();
        }
        if (GUILayout.Button("Save"))
        {
            SaveSOs();
        }
        // Reset button to reset startingPoly
        if (GUILayout.Button("Reset"))
        {
            poly = new Polyhedron(startingMesh);//  startingPoly = Polyhedron.ComputeIcosahedron();
            RecomputePolyMesh();
            
        }
        GUI.enabled = isProcessing;
        if (GUILayout.Button("Cancel"))//,buttonStyle))
        {
            CancelTask();
        }
        GUI.enabled = true;

        if (isProcessing)
            DrawOscillatingBox(new Vector2(500, 500));


        void DrawOscillatingBox(Vector2 pos)
        {
            // Update box size using a sine wave
            float frac = Mathf.Sin(Time.time * speed) * 0.4f + 0.6f; // Oscillates between 0.2 and 1

            // Calculate dimensions
            float width = boxSize * (1f+frac);
            float height = boxSize * 0.5f;// * (1f - frac); // Opposite phase

            // Draw the box centered on the screen
            Rect boxRect = new Rect(pos.x-(width/2f),pos.y, width, height);
            GUI.Box(boxRect, "Processing...");
        }
    }
    public float boxSize = 200f; // Box size multiplier (oscillates)
    public float speed = 1f;     // Oscillation speed
}
