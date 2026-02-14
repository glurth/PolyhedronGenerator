using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Cysharp.Threading.Tasks;
using EyE.Threading;
using System.Text.RegularExpressions;
public class CommandLineProcessor
{
    private readonly Dictionary<string, Func<UniTask>> operations;
    private readonly CancelBoolRef cancelRef;

    public CommandLineProcessor(
        Dictionary<string, Func<UniTask>> operations,
        CancelBoolRef getCancelRef)
    {
        this.operations = operations;
        this.cancelRef = getCancelRef;
    }

    public async UniTask ExecuteSequence(string commandInput)
    {
        List<(string, int)> operationsToRun = ParseCommands(commandInput);
        int opCount = 0;
        foreach ((string cmd, int count) in operationsToRun)
        {
            if (operations.TryGetValue(cmd, out Func<UniTask> op))
            {
                for (int i = 0; i < count; i++)
                {
                    Debug.Log("["+opCount+"] Executing Command: "+cmd);
                    await op();
                }
                opCount++;
            }
            else
            {
                Debug.LogWarning($"Unknown command: {cmd}");
            }

        }
    }

    private List<(string, int)> ParseCommands(string input)
    {
        List<(string, int)> result = new List<(string, int)>();
        MatchCollection matches = Regex.Matches(input, @"([a-zA-Z]+)\s*(\d*)");
        Debug.Log("Commands detected: " + matches.Count);
        foreach (Match match in matches)
        {
            string cmd = match.Groups[1].Value.ToLower();
            int count = 1;
            if (!string.IsNullOrEmpty(match.Groups[2].Value)) 
                count = int.Parse(match.Groups[2].Value);
            string key = ResolveCommandKey(cmd);
            result.Add((key, count));
        }
        return result;
    }

    private string ResolveCommandKey(string input)
    {
        foreach (string key in operations.Keys)
        {
            if (key.StartsWith(input)) return key;
        }
        Debug.Log("unable to resolve command for: " + input);
        return input;
    }
}

